#!/usr/bin/env python3
"""
train_rs25.py

Train the joint detect+type base-pair model on a single graph set (default the
RS25 graphs) with a structure-level **train / val / test** split. The epoch is
selected by validation detection-F1 and final metrics are reported on the
held-out **test** split (which the model never trained on or selected against).

Usage:
    python train_rs25.py [graphs_dir] --epochs 30
"""
from __future__ import annotations

import argparse
import copy
import sys
from pathlib import Path

import numpy as np
import torch

import train_pair_detect_type as dt
from rna_pairs import per_class_prf


def split_indices(n: int, val_frac: float, test_frac: float, seed: int):
    """Structure-level split -> (train_idx, val_idx, test_idx)."""
    perm = torch.randperm(n, generator=torch.Generator().manual_seed(seed)).tolist()
    n_test = int(round(test_frac * n))
    n_val = int(round(val_frac * n))
    return perm[n_test + n_val:], perm[n_test:n_test + n_val], perm[:n_test]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("graphs", nargs="?", default=Path("rs25/graphs"), type=Path)
    p.add_argument("--epochs", type=int, default=30)
    p.add_argument("--val-frac", type=float, default=0.15)
    p.add_argument("--test-frac", type=float, default=0.15)
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--lr", type=float, default=5e-3)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    torch.manual_seed(args.seed)
    gen = torch.Generator().manual_seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    data, class_names, _ = dt.load_dataset(args.graphs, seed=args.seed)
    tr_i, va_i, te_i = split_indices(len(data), args.val_frac, args.test_frac, args.seed)
    train = [data[i] for i in tr_i]
    val = [data[i] for i in va_i]
    test = [data[i] for i in te_i]
    print(f"[rs25] {len(data)} graphs -> train {len(train)} / val {len(val)} / test {len(test)} "
          f"| {len(class_names)} classes", file=sys.stderr)

    model = dt.PairTypeGNN(num_classes=len(class_names)).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)

    best_f1, best_state = -1.0, None
    for epoch in range(1, args.epochs + 1):
        dt.run_epoch(model, train, train=True, optimizer=opt, batch_size=args.batch_size, gen=gen, device=device)
        vm = dt.run_epoch(model, val, train=False, optimizer=opt, batch_size=args.batch_size, gen=gen, device=device)
        if vm["det_F1"] > best_f1:                      # model selection on validation
            best_f1, best_state = vm["det_F1"], copy.deepcopy(model.state_dict())
        print(f"epoch {epoch:3d} | val detF1 {vm['det_F1']:.3f} exact {vm['exact_acc']:.3f} "
              f"typing {vm['type_acc_given_detected']:.3f}")

    model.load_state_dict(best_state)                   # restore best-val checkpoint
    conf = dt.confusion(model, test)
    tm = dt.metrics_from_confusion(conf)
    print(f"\n[rs25] TEST (checkpoint = best val detF1 {best_f1:.3f}):", file=sys.stderr)
    print(f"  detection F1 {tm['det_F1']:.3f} | exact (detect+type) {tm['exact_acc']:.3f} "
          f"| typing | detected {tm['type_acc_given_detected']:.3f}")
    prec, rec, f1, sup = per_class_prf(conf)
    order = [i for i in np.argsort(-sup.numpy()) if sup[i] > 0]
    print("  per-class F1 (test):",
          ", ".join(f"{class_names[i]}={float(f1[i]):.2f}(n={int(sup[i])})" for i in order))
    return 0


if __name__ == "__main__":
    sys.exit(main())
