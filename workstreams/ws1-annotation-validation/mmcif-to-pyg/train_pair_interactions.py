#!/usr/bin/env python3
"""
train_pair_interactions.py

Predict the **interaction type of base pairs** (edge classification).

For every annotated base pair we predict its Leontis-Westhof family (cWW, tWS,
cSH, ...). To keep it well-posed (no label leakage), the message-passing graph
tells the model *which* residues pair via a single generic ``PAIR`` relation and
the backbone (``B53``/``B35``) -- but NOT the specific family. From node
identities + this topology the RGCN learns the geometry class of each pair.

  relations (message passing): B53, B35, PAIR        (3)
  node features              : one-hot nt {A,C,G,U,N}
  per-pair target            : LW family of the lo->hi edge (lo<hi by node index)

Usage:
    python train_pair_interactions.py graphs/ --epochs 25
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import torch
from torch_geometric.data import Batch

from rna_pairs import (NT, NT_IDX, PairData, PairGNN, classification_metrics, is_canonical,
                       load_graphs, parse_graph, train_eval_classification)

PAIR_REL = 2          # message-passing relation id for "these residues pair"
NUM_RELATIONS = 3     # B53, B35, PAIR


def load_dataset(graph_dir: Path):
    """rnaglib JSON graphs -> (list[PairData], lw_class_vocab)."""
    _, graphs = load_graphs(graph_dir)
    parsed = [parse_graph(G) for G in graphs]

    fam_counts: Counter = Counter()
    for pg in parsed:
        for _, _, fam in pg.base_pairs:
            fam_counts[fam] += 1
    lw_vocab = {fam: i for i, fam in enumerate(sorted(fam_counts))}

    data_list: list[PairData] = []
    for pg in parsed:
        if not pg.base_pairs:
            continue
        src, dst, rel = list(pg.bb_src), list(pg.bb_dst), list(pg.bb_rel)
        lo_idx, hi_idx, py, nt = [], [], [], []
        for lo, hi, fam in pg.base_pairs:
            src += [lo, hi]; dst += [hi, lo]; rel += [PAIR_REL, PAIR_REL]   # pair edge, both directions
            lo_idx.append(lo); hi_idx.append(hi); py.append(lw_vocab[fam])
            nt.append([int(pg.x[lo].argmax()), int(pg.x[hi].argmax())])     # base identities, for canonicity

        data = PairData(
            x=pg.x,
            edge_index=torch.tensor([src, dst], dtype=torch.long),
            edge_type=torch.tensor(rel, dtype=torch.long),
            pair_index=torch.tensor([lo_idx, hi_idx], dtype=torch.long),
            pair_y=torch.tensor(py, dtype=torch.long),
            pair_nt=torch.tensor(nt, dtype=torch.long),
        )
        data.num_nodes = pg.num_nodes
        data_list.append(data)
    return data_list, lw_vocab


class PairRGCN(PairGNN):
    def __init__(self, num_classes: int, hidden: int = 64):
        super().__init__(out_dim=num_classes, num_relations=NUM_RELATIONS, hidden=hidden, use_geom=False)


def _noncanonical_recall(model, data_list, device, include_wobble: bool = False) -> float:
    """Base-aware: among truly non-canonical pairs (cWW + A-U/G-C excluded), the
    fraction whose LW family the model predicts correctly. Needs model._inv_vocab."""
    inv = getattr(model, "_inv_vocab", None)
    if inv is None:
        return 0.0
    model.eval()
    correct = total = 0
    with torch.no_grad():
        for d in data_list:
            b = Batch.from_data_list([d]).to(device)
            pred = model(b.x, b.edge_index, b.edge_type, b.pair_index).argmax(1).cpu()
            for k in range(d.pair_y.numel()):
                fam = inv[int(d.pair_y[k])]
                b1, b2 = NT[int(d.pair_nt[k, 0])], NT[int(d.pair_nt[k, 1])]
                if not is_canonical(fam, b1, b2, include_wobble):
                    total += 1
                    correct += int(pred[k] == d.pair_y[k])
    return correct / max(total, 1)


def run_epoch(model, data_list, *, train, optimizer, batch_size, gen, device) -> dict:
    conf, loss = train_eval_classification(model, data_list, train=train, optimizer=optimizer,
                                           batch_size=batch_size, gen=gen, device=device)
    metrics = classification_metrics(conf)
    metrics["loss"] = loss
    # base-aware recall over genuinely non-canonical pairs (val only; it costs an extra pass)
    metrics["noncanonical_recall"] = 0.0 if train else _noncanonical_recall(model, data_list, device)
    return metrics


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("graphs", type=Path)
    p.add_argument("--epochs", type=int, default=25)
    p.add_argument("--batch-size", type=int, default=32)
    p.add_argument("--hidden", type=int, default=64)
    p.add_argument("--lr", type=float, default=1e-2)
    p.add_argument("--val-frac", type=float, default=0.2)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    torch.manual_seed(args.seed)
    gen = torch.Generator().manual_seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    data_list, lw_vocab = load_dataset(args.graphs)
    inv = {v: k for k, v in lw_vocab.items()}
    print(f"[pairs] {len(data_list)} graphs | LW classes: {lw_vocab}", file=sys.stderr)

    perm = torch.randperm(len(data_list), generator=gen).tolist()
    n_val = max(1, int(len(data_list) * args.val_frac))
    val = [data_list[i] for i in perm[:n_val]]
    train = [data_list[i] for i in perm[n_val:]] or val

    counts = Counter()
    for d in train:
        counts.update(d.pair_y.tolist())
    baseline = max(counts.values()) / sum(counts.values())
    print(f"[pairs] class counts {{{', '.join(f'{inv[c]}:{n}' for c, n in counts.most_common())}}}", file=sys.stderr)
    print(f"[pairs] majority baseline ~{baseline:.3f}", file=sys.stderr)

    model = PairRGCN(num_classes=len(lw_vocab), hidden=args.hidden).to(device)
    model._inv_vocab = inv          # enables base-aware non-canonical recall
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)

    for epoch in range(1, args.epochs + 1):
        run_epoch(model, train, train=True, optimizer=opt, batch_size=args.batch_size, gen=gen, device=device)
        va = run_epoch(model, val, train=False, optimizer=opt, batch_size=args.batch_size, gen=gen, device=device)
        print(f"epoch {epoch:3d} | val micro {va['micro_acc']:.3f} | balanced {va['balanced_acc']:.3f} "
              f"| macro-F1 {va['macro_f1']:.3f} | non-canon recall {va['noncanonical_recall']:.3f}")
    print(f"[pairs] final | micro {va['micro_acc']:.3f} (majority baseline {baseline:.3f}) "
          f"| balanced {va['balanced_acc']:.3f} | macro-F1 {va['macro_f1']:.3f} "
          f"| non-canon recall {va['noncanonical_recall']:.3f}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
