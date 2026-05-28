#!/usr/bin/env python3
"""
train_dnatco.py

Train a simple model on the rnaglib-format graphs produced by dnatco_to_graph.py.

Task: **masked nucleotide prediction** (self-supervised). We hide the identity of
a fraction of residues (replace with a MASK token) and ask a Relational GCN to
recover A/C/G/U from the graph context. Because base-paired partners are
complementary (cWW G-C, A-U) and backbone neighbours are informative, the model
can do meaningfully better than the majority-class baseline.

  - node input : one-hot over {A, C, G, U, MASK}
  - edge types : Leontis-Westhof families + backbone (B53/B35), used as RGCN relations
  - target     : true base at masked, non-modified residues

Usage:
    python train_dnatco.py graphs/ --epochs 20
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import torch
import torch.nn.functional as F
from torch_geometric.data import Batch, Data
from torch_geometric.nn import RGCNConv

NT = ["A", "C", "G", "U"]
NT_IDX = {c: i for i, c in enumerate(NT)}
MASK = len(NT)          # input token index for "hidden"
IN_DIM = len(NT) + 1    # A,C,G,U,MASK


def load_dataset(graph_dir: Path):
    """Read rnaglib JSON graphs -> (list[Data], edge_vocab)."""
    from rnaglib.utils import load_graph

    files = sorted(graph_dir.glob("*.json"))
    if not files:
        raise SystemExit(f"error: no .json graphs in {graph_dir}")

    raw = []
    edge_types: set[str] = set()
    for f in files:
        G = load_graph(str(f))
        edge_types.update(d["LW"] for *_, d in G.edges(data=True))
        raw.append(G)
    edge_vocab = {lw: i for i, lw in enumerate(sorted(edge_types))}

    data_list = []
    for G in raw:
        nodes = sorted(G.nodes, key=lambda n: (G.nodes[n]["chain_id"], G.nodes[n]["index"]))
        idx = {n: i for i, n in enumerate(nodes)}
        y = torch.tensor([NT_IDX.get(G.nodes[n]["nt_code"], -1) for n in nodes])
        src, dst, et = [], [], []
        for u, v, d in G.edges(data=True):
            src.append(idx[u]); dst.append(idx[v]); et.append(edge_vocab[d["LW"]])
        data = Data(
            edge_index=torch.tensor([src, dst], dtype=torch.long) if src else torch.empty((2, 0), dtype=torch.long),
            edge_type=torch.tensor(et, dtype=torch.long),
            y=y,
        )
        data.num_nodes = len(nodes)
        if (y >= 0).any():           # keep graphs with at least one standard base
            data_list.append(data)
    return data_list, edge_vocab


def mask_inputs(data: Data, mask_rate: float, generator: torch.Generator):
    """Return a masked copy of `data` (x one-hot with MASK) + which nodes are scored."""
    n = data.num_nodes
    valid = data.y >= 0
    rand = torch.rand(n, generator=generator) < mask_rate
    masked = valid & rand
    revealed = valid & ~masked

    x = torch.zeros(n, IN_DIM)
    x[revealed, data.y[revealed]] = 1.0      # show true base
    x[~revealed, MASK] = 1.0                 # masked valid + modified/unknown -> MASK

    out = Data(x=x, edge_index=data.edge_index, edge_type=data.edge_type, y=data.y)
    out.num_nodes = n
    out.score_mask = masked
    return out


class RGCN(torch.nn.Module):
    def __init__(self, num_relations: int, hidden: int = 64):
        super().__init__()
        self.conv1 = RGCNConv(IN_DIM, hidden, num_relations)
        self.conv2 = RGCNConv(hidden, hidden, num_relations)
        self.lin = torch.nn.Linear(hidden, len(NT))

    def forward(self, x, edge_index, edge_type):
        x = F.relu(self.conv1(x, edge_index, edge_type))
        x = F.relu(self.conv2(x, edge_index, edge_type))
        return self.lin(x)


def run_epoch(model, data_list, edge_vocab, *, train, optimizer, mask_rate, batch_size, gen, device):
    model.train(train)
    order = torch.randperm(len(data_list), generator=gen).tolist() if train else range(len(data_list))
    total_loss = correct = scored = 0
    for start in range(0, len(data_list), batch_size):
        chunk = [mask_inputs(data_list[i], mask_rate, gen) for i in (order[start:start + batch_size] if train else list(order)[start:start + batch_size])]
        batch = Batch.from_data_list(chunk).to(device)
        if batch.score_mask.sum() == 0:
            continue
        with torch.set_grad_enabled(train):
            logits = model(batch.x, batch.edge_index, batch.edge_type)
            m = batch.score_mask
            loss = F.cross_entropy(logits[m], batch.y[m])
            if train:
                optimizer.zero_grad(); loss.backward(); optimizer.step()
        total_loss += float(loss.detach()) * int(m.sum())
        correct += int((logits[m].argmax(1) == batch.y[m]).sum())
        scored += int(m.sum())
    return total_loss / max(scored, 1), correct / max(scored, 1)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("graphs", type=Path, help="Directory of .json graphs.")
    p.add_argument("--epochs", type=int, default=20)
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--mask-rate", type=float, default=0.15)
    p.add_argument("--hidden", type=int, default=64)
    p.add_argument("--lr", type=float, default=1e-2)
    p.add_argument("--val-frac", type=float, default=0.2)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    torch.manual_seed(args.seed)
    gen = torch.Generator().manual_seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    data_list, edge_vocab = load_dataset(args.graphs)
    print(f"[train] {len(data_list)} graphs | edge relations: {edge_vocab}", file=sys.stderr)

    n_val = max(1, int(len(data_list) * args.val_frac))
    perm = torch.randperm(len(data_list), generator=gen).tolist()
    val = [data_list[i] for i in perm[:n_val]]
    train = [data_list[i] for i in perm[n_val:]] or val

    # majority-class baseline on val (most frequent base among scored positions)
    counts = torch.zeros(len(NT))
    for d in train:
        for c in d.y[d.y >= 0]:
            counts[c] += 1
    baseline = float(counts.max() / counts.sum())
    print(f"[train] base frequencies {counts.tolist()} | majority baseline ~{baseline:.3f}", file=sys.stderr)

    model = RGCN(num_relations=len(edge_vocab), hidden=args.hidden).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)

    for epoch in range(1, args.epochs + 1):
        tr_loss, tr_acc = run_epoch(model, train, edge_vocab, train=True, optimizer=opt,
                                    mask_rate=args.mask_rate, batch_size=args.batch_size, gen=gen, device=device)
        va_loss, va_acc = run_epoch(model, val, edge_vocab, train=False, optimizer=opt,
                                    mask_rate=args.mask_rate, batch_size=args.batch_size, gen=gen, device=device)
        print(f"epoch {epoch:3d} | train loss {tr_loss:.3f} acc {tr_acc:.3f} "
              f"| val loss {va_loss:.3f} acc {va_acc:.3f}")
    print(f"[train] final val accuracy {va_acc:.3f} (baseline {baseline:.3f})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
