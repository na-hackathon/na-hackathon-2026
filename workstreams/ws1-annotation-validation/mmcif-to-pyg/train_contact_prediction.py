#!/usr/bin/env python3
"""
train_contact_prediction.py

Predict **which residues base-pair** (contact / link prediction) using the
experimental geometry stored in the graphs (xyz_C1p / xyz_glyN / xyz_P).

This is the realistic structural task: from sequence + backbone + 3D geometry,
decide whether two nucleotides form a base pair.

Design (kept honest):
  - message passing uses the **backbone only** (B53/B35) -- the base-pair edges
    are the labels, so they are NOT given to the network.
  - geometry enters as rotation/translation-invariant **distance features** per
    candidate pair: C1'-C1', glyN-glyN, P-P, plus sequence separation.
  - **hard negatives**: non-paired pairs that are spatially close (within a
    distance cutoff, via a KD-tree) -- so the task can't be solved by distance
    alone, the model must use identity + orientation + context.

Reports precision / recall / F1 / accuracy against a distance-only baseline.

Usage:
    python train_contact_prediction.py graphs/ --epochs 20 [--no-geom]
(`graphs/` must contain coordinates -- build with dnatco_to_graph.py, coords on by default.)
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import torch

from rna_pairs import (PairData, PairGNN, binary_prf, geom_features, load_graphs,
                       parse_graph, sample_hard_negatives, train_eval_binary)

NUM_RELATIONS = 2     # backbone only (B53, B35)


def load_dataset(graph_dir: Path, neg_ratio: float = 1.0, cutoff: float = 18.0, seed: int = 0):
    """rnaglib JSON graphs (with coords) -> list[PairData] for contact prediction."""
    _, graphs = load_graphs(graph_dir)
    rng = np.random.default_rng(seed)

    data_list, skipped = [], 0
    for G in graphs:
        pg = parse_graph(G)
        if sum(c is not None for c in pg.c1) < 4:
            skipped += 1; continue
        pos = {(lo, hi) for lo, hi, _ in pg.base_pairs
               if pg.c1[lo] is not None and pg.c1[hi] is not None}
        if not pos:
            skipped += 1; continue

        negs = sample_hard_negatives(pg.c1, pg.seq, pg.chain, pos, cutoff,
                                     int(round(neg_ratio * len(pos))), rng)
        pairs = list(pos) + negs
        labels = [1.0] * len(pos) + [0.0] * len(negs)
        feats, dc1 = [], []
        for a, b in pairs:
            f, d = geom_features(a, b, pg.seq, pg.chain, pg.c1, pg.gly, pg.pp)
            feats.append(f); dc1.append(d)

        data = PairData(
            x=pg.x,
            edge_index=torch.tensor([pg.bb_src, pg.bb_dst], dtype=torch.long) if pg.bb_src else torch.empty((2, 0), dtype=torch.long),
            edge_type=torch.tensor(pg.bb_rel, dtype=torch.long),
            pair_index=torch.tensor([[a for a, _ in pairs], [b for _, b in pairs]], dtype=torch.long),
            pair_geom=torch.tensor(feats, dtype=torch.float),
            pair_y=torch.tensor(labels, dtype=torch.float),
            pair_dc1=torch.tensor(dc1, dtype=torch.float),
        )
        data.num_nodes = pg.num_nodes
        data_list.append(data)
    if skipped:
        print(f"[contact] skipped {skipped} graphs (no coords or no pairs)", file=sys.stderr)
    return data_list


class ContactGNN(PairGNN):
    def __init__(self, hidden: int = 64, use_geom: bool = True):
        super().__init__(out_dim=1, num_relations=NUM_RELATIONS, hidden=hidden, use_geom=use_geom)


def run_epoch(model, data_list, *, train, optimizer, batch_size, gen, device) -> dict:
    pred, y, loss = train_eval_binary(model, data_list, train=train, optimizer=optimizer,
                                      batch_size=batch_size, gen=gen, device=device)
    metrics = binary_prf(pred, y)
    metrics["loss"] = loss
    return metrics


def distance_baseline(data_list, thresh: float = 12.5) -> dict:
    """Predict pair iff C1'-C1' distance <= thresh (Angstrom)."""
    preds = [(d.pair_dc1 <= thresh).long() for d in data_list]
    ys = [d.pair_y.long() for d in data_list]
    return binary_prf(torch.cat(preds), torch.cat(ys))


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("graphs", type=Path)
    p.add_argument("--epochs", type=int, default=20)
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--hidden", type=int, default=64)
    p.add_argument("--lr", type=float, default=5e-3)
    p.add_argument("--neg-ratio", type=float, default=1.0)
    p.add_argument("--cutoff", type=float, default=18.0, help="KD-tree distance (A) for hard negatives.")
    p.add_argument("--no-geom", action="store_true",
                   help="Ablation: drop geometric distance features (sequence + backbone only).")
    p.add_argument("--val-frac", type=float, default=0.2)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    torch.manual_seed(args.seed)
    gen = torch.Generator().manual_seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    data_list = load_dataset(args.graphs, neg_ratio=args.neg_ratio, cutoff=args.cutoff, seed=args.seed)
    print(f"[contact] {len(data_list)} graphs | "
          f"{sum(int(d.pair_y.sum()) for d in data_list)} positive pairs", file=sys.stderr)

    perm = torch.randperm(len(data_list), generator=gen).tolist()
    n_val = max(1, int(len(data_list) * args.val_frac))
    val = [data_list[i] for i in perm[:n_val]]
    train = [data_list[i] for i in perm[n_val:]] or val

    base = distance_baseline(val)
    print(f"[contact] distance-only baseline (val): "
          f"F1 {base['f1']:.3f} (P {base['precision']:.3f} R {base['recall']:.3f})", file=sys.stderr)

    model = ContactGNN(hidden=args.hidden, use_geom=not args.no_geom).to(device)
    mode = "sequence+backbone only (no geom)" if args.no_geom else "with geometry"
    print(f"[contact] model: {mode}", file=sys.stderr)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    for epoch in range(1, args.epochs + 1):
        run_epoch(model, train, train=True, optimizer=opt, batch_size=args.batch_size, gen=gen, device=device)
        m = run_epoch(model, val, train=False, optimizer=opt, batch_size=args.batch_size, gen=gen, device=device)
        print(f"epoch {epoch:3d} | val acc {m['acc']:.3f} | P {m['precision']:.3f} "
              f"R {m['recall']:.3f} F1 {m['f1']:.3f}")
    print(f"[contact] final val F1 {m['f1']:.3f} ({mode}) vs distance-only baseline F1 {base['f1']:.3f}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
