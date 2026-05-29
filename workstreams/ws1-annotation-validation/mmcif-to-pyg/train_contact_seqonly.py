#!/usr/bin/env python3
"""
train_contact_seqonly.py

Coordinate-FREE base-pair contact prediction: decide whether two nucleotides
base-pair from **sequence + backbone topology only** -- no coordinates anywhere,
*including candidate generation*. This is the hard structure-prediction setting,
unlike train_contact_prediction.py, which is handed the 3D geometry.

Candidates per structure:
  - positives = all annotated base pairs
  - negatives = randomly sampled non-pairs with sequence separation >= --min-sep
    (a sequence-based prior, NOT a spatial/KD-tree filter)

Model: RGCN over the backbone (B53/B35) with one-hot nucleotide node features;
the pair head sees only [h_i, h_j, h_i*h_j].

`--with-geom` adds the C1'/glyN/P distance features on the SAME sequence-sampled
candidate set (requires coords) -- a fair "does geometry help" ablation that holds
the candidate set fixed (so the only change is the model's features, not the negatives).

Usage:
    python train_contact_seqonly.py graphs/ --epochs 25
    python train_contact_seqonly.py graphs/ --epochs 25 --with-geom
"""
from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import numpy as np
import torch

from rna_pairs import (GEOM_DIM, PairData, PairGNN, binary_prf, geom_features,
                       load_graphs, parse_graph, train_eval_binary)

NUM_RELATIONS = 2     # backbone only (B53, B35)


def sample_seq_negatives(num_nodes, seq, chain, positives, n_neg, min_sep, rng):
    """Random non-pair (lo, hi) with sequence separation >= min_sep (no geometry)."""
    negs: set = set()
    attempts, cap = 0, max(2000, n_neg * 60)
    while len(negs) < n_neg and attempts < cap:
        attempts += 1
        i, j = int(rng.integers(num_nodes)), int(rng.integers(num_nodes))
        if i == j:
            continue
        lo, hi = (i, j) if i < j else (j, i)
        if (lo, hi) in positives or (lo, hi) in negs:
            continue
        if chain[lo] == chain[hi] and abs(seq[lo] - seq[hi]) < min_sep:
            continue
        negs.add((lo, hi))
    return list(negs)


def load_dataset(graph_dir: Path, neg_ratio: float = 1.0, min_sep: int = 3,
                 seed: int = 0, with_geom: bool = False):
    """rnaglib JSON graphs -> list[PairData] with sequence-sampled candidates."""
    _, graphs = load_graphs(graph_dir)
    rng = np.random.default_rng(seed)

    data_list, skipped = [], 0
    for G in graphs:
        pg = parse_graph(G)
        pos = {(lo, hi) for lo, hi, _ in pg.base_pairs}
        if not pos or pg.num_nodes < 6:
            skipped += 1; continue
        if with_geom and sum(c is not None for c in pg.c1) < 4:
            skipped += 1; continue

        negs = sample_seq_negatives(pg.num_nodes, pg.seq, pg.chain, pos,
                                    int(round(neg_ratio * len(pos))), min_sep, rng)
        pairs = list(pos) + negs
        labels = [1.0] * len(pos) + [0.0] * len(negs)

        kw = {}
        if with_geom:
            feats = []
            for a, b in pairs:
                if pg.c1[a] is not None and pg.c1[b] is not None:
                    feats.append(geom_features(a, b, pg.seq, pg.chain, pg.c1, pg.gly, pg.pp)[0])
                else:
                    feats.append([0.0] * GEOM_DIM)
            kw["pair_geom"] = torch.tensor(feats, dtype=torch.float)

        data = PairData(
            x=pg.x,
            edge_index=torch.tensor([pg.bb_src, pg.bb_dst], dtype=torch.long) if pg.bb_src else torch.empty((2, 0), dtype=torch.long),
            edge_type=torch.tensor(pg.bb_rel, dtype=torch.long),
            pair_index=torch.tensor([[a for a, _ in pairs], [b for _, b in pairs]], dtype=torch.long),
            pair_y=torch.tensor(labels, dtype=torch.float),
            **kw,
        )
        data.num_nodes = pg.num_nodes
        data_list.append(data)
    if skipped:
        print(f"[seqonly] skipped {skipped} graphs (no pairs / too small)", file=sys.stderr)
    return data_list


class SeqContactGNN(PairGNN):
    def __init__(self, hidden: int = 64, use_geom: bool = False):
        super().__init__(out_dim=1, num_relations=NUM_RELATIONS, hidden=hidden, use_geom=use_geom)


def run_epoch(model, data_list, *, train, optimizer, batch_size, gen, device) -> dict:
    pred, y, loss = train_eval_binary(model, data_list, train=train, optimizer=optimizer,
                                      batch_size=batch_size, gen=gen, device=device)
    metrics = binary_prf(pred, y)
    metrics["loss"] = loss
    return metrics


def split_indices(n, val_frac, test_frac, seed):
    perm = torch.randperm(n, generator=torch.Generator().manual_seed(seed)).tolist()
    n_test = int(round(test_frac * n))
    n_val = int(round(val_frac * n))
    return perm[n_test + n_val:], perm[n_test:n_test + n_val], perm[:n_test]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("graphs", type=Path)
    p.add_argument("--epochs", type=int, default=25)
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--lr", type=float, default=5e-3)
    p.add_argument("--neg-ratio", type=float, default=1.0)
    p.add_argument("--min-sep", type=int, default=3, help="min sequence separation for negatives.")
    p.add_argument("--with-geom", action="store_true", help="add geometry features on the same candidates.")
    p.add_argument("--val-frac", type=float, default=0.15)
    p.add_argument("--test-frac", type=float, default=0.15)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    torch.manual_seed(args.seed)
    gen = torch.Generator().manual_seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    data = load_dataset(args.graphs, neg_ratio=args.neg_ratio, min_sep=args.min_sep,
                        seed=args.seed, with_geom=args.with_geom)
    tr_i, va_i, te_i = split_indices(len(data), args.val_frac, args.test_frac, args.seed)
    train = [data[i] for i in tr_i]; val = [data[i] for i in va_i]; test = [data[i] for i in te_i]
    mode = "sequence + geometry (same candidates)" if args.with_geom else "sequence only (no coordinates)"
    print(f"[seqonly] {len(data)} graphs -> train {len(train)} / val {len(val)} / test {len(test)} | mode: {mode}", file=sys.stderr)

    model = SeqContactGNN(use_geom=args.with_geom).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    best_f1, best_state = -1.0, None
    for epoch in range(1, args.epochs + 1):
        run_epoch(model, train, train=True, optimizer=opt, batch_size=args.batch_size, gen=gen, device=device)
        vm = run_epoch(model, val, train=False, optimizer=opt, batch_size=args.batch_size, gen=gen, device=device)
        if vm["f1"] > best_f1:
            best_f1, best_state = vm["f1"], copy.deepcopy(model.state_dict())
        print(f"epoch {epoch:3d} | val F1 {vm['f1']:.3f} (P {vm['precision']:.3f} R {vm['recall']:.3f})")
    model.load_state_dict(best_state)
    tm = run_epoch(model, test, train=False, optimizer=opt, batch_size=args.batch_size, gen=gen, device=device)
    print(f"\n[seqonly] TEST ({mode}) | F1 {tm['f1']:.3f} | P {tm['precision']:.3f} | R {tm['recall']:.3f} "
          f"| acc {tm['acc']:.3f}  (best val F1 {best_f1:.3f})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
