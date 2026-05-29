#!/usr/bin/env python3
"""
train_pair_detect_type.py

Joint **detection + typing** of base pairs: for each candidate residue pair,
predict either NONE (not a base pair) or its Leontis-Westhof family.

This unifies the two earlier models:
  - Model 2 (contact prediction) only said pair / no-pair.
  - Model 1 (interaction type) only typed pairs it was *told* were pairs.
Here the model is NOT told which residues pair; from sequence + backbone +
geometry it must both find the pairs AND say what type each is.

Classes: {NONE=0} + the LW families present (cWW, tWS, ...).
Same graph inputs as train_contact_prediction (backbone-only message passing,
geometric distance features, hard spatially-close negatives).

Reports:
  - detection P / R / F1       (NONE vs any pair)
  - exact accuracy on true pairs (detected AND typed correctly)
  - typing accuracy | detected (of pairs it found, fraction typed correctly)

Usage:
    python train_pair_detect_type.py graphs/ --epochs 20 [--no-geom]
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import torch

# helpers re-exported below are the public API the notebook imports as `dt.*`
from rna_pairs import (  # noqa: F401
    PairData, PairGNN, confusion, geom_features, graph_node_order, load_graphs,
    parse_graph, per_class_prf, per_structure_metrics, predict_one,
    sample_hard_negatives, train_eval_classification,
    detection_typing_metrics as metrics_from_confusion,
)

NUM_RELATIONS = 2     # backbone only (B53, B35)


def load_dataset(graph_dir: Path, neg_ratio: float = 1.0, cutoff: float = 18.0,
                 seed: int = 0, class_names: list | None = None):
    """rnaglib JSON graphs (with coords) -> (list[PairData], class_names, kept_files).

    pair_y: 0 = NONE (not a pair), k>0 = LW family index. Negatives are
    spatially-close non-pairs (hard negatives), same as contact prediction.

    Pass `class_names` (["NONE", fam, ...]) to encode labels against a FIXED
    vocabulary -- required to evaluate a trained model on a different dataset so
    class indices line up. Pairs whose family is not in the vocab are skipped.
    """
    files, graphs = load_graphs(graph_dir)
    parsed = [parse_graph(G) for G in graphs]
    rng = np.random.default_rng(seed)

    if class_names is None:
        fam_counts: Counter = Counter()
        for pg in parsed:
            for _, _, fam in pg.base_pairs:
                fam_counts[fam] += 1
        class_names = ["NONE"] + sorted(fam_counts)
    lw_class = {fam: i + 1 for i, fam in enumerate(class_names[1:])}

    data_list, kept_files, skipped = [], [], 0
    for fpath, pg in zip(files, parsed):
        if sum(c is not None for c in pg.c1) < 4:
            skipped += 1; continue

        pos = {}
        for lo, hi, fam in pg.base_pairs:
            if pg.c1[lo] is None or pg.c1[hi] is None or fam not in lw_class:
                continue
            pos[(lo, hi)] = lw_class[fam]
        if not pos:
            skipped += 1; continue

        negs = sample_hard_negatives(pg.c1, pg.seq, pg.chain, set(pos),
                                     cutoff, int(round(neg_ratio * len(pos))), rng)
        pairs = list(pos) + negs
        labels = list(pos.values()) + [0] * len(negs)
        feats = [geom_features(a, b, pg.seq, pg.chain, pg.c1, pg.gly, pg.pp)[0] for a, b in pairs]

        data = PairData(
            x=pg.x,
            edge_index=torch.tensor([pg.bb_src, pg.bb_dst], dtype=torch.long) if pg.bb_src else torch.empty((2, 0), dtype=torch.long),
            edge_type=torch.tensor(pg.bb_rel, dtype=torch.long),
            pair_index=torch.tensor([[a for a, _ in pairs], [b for _, b in pairs]], dtype=torch.long),
            pair_geom=torch.tensor(feats, dtype=torch.float),
            pair_y=torch.tensor(labels, dtype=torch.long),
        )
        data.num_nodes = pg.num_nodes
        data_list.append(data); kept_files.append(fpath)
    if skipped:
        print(f"[detect+type] skipped {skipped} graphs (no coords or no pairs)", file=sys.stderr)
    return data_list, class_names, kept_files


class PairTypeGNN(PairGNN):
    def __init__(self, num_classes: int, hidden: int = 64, use_geom: bool = True):
        super().__init__(out_dim=num_classes, num_relations=NUM_RELATIONS, hidden=hidden, use_geom=use_geom)


def run_epoch(model, data_list, *, train, optimizer, batch_size, gen, device) -> dict:
    conf, loss = train_eval_classification(model, data_list, train=train, optimizer=optimizer,
                                           batch_size=batch_size, gen=gen, device=device)
    metrics = metrics_from_confusion(conf)
    metrics["loss"] = loss
    return metrics


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("graphs", type=Path)
    p.add_argument("--epochs", type=int, default=20)
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--hidden", type=int, default=64)
    p.add_argument("--lr", type=float, default=5e-3)
    p.add_argument("--neg-ratio", type=float, default=1.0)
    p.add_argument("--cutoff", type=float, default=18.0)
    p.add_argument("--no-geom", action="store_true")
    p.add_argument("--val-frac", type=float, default=0.2)
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    torch.manual_seed(args.seed)
    gen = torch.Generator().manual_seed(args.seed)
    device = "cuda" if torch.cuda.is_available() else "cpu"

    data_list, class_names, _ = load_dataset(args.graphs, neg_ratio=args.neg_ratio, cutoff=args.cutoff, seed=args.seed)
    print(f"[detect+type] {len(data_list)} graphs | {len(class_names)} classes "
          f"(NONE + {len(class_names) - 1} LW families)", file=sys.stderr)

    perm = torch.randperm(len(data_list), generator=gen).tolist()
    n_val = max(1, int(len(data_list) * args.val_frac))
    val = [data_list[i] for i in perm[:n_val]]
    train = [data_list[i] for i in perm[n_val:]] or val

    model = PairTypeGNN(num_classes=len(class_names), use_geom=not args.no_geom).to(device)
    opt = torch.optim.Adam(model.parameters(), lr=args.lr)
    for epoch in range(1, args.epochs + 1):
        run_epoch(model, train, train=True, optimizer=opt, batch_size=args.batch_size, gen=gen, device=device)
        m = run_epoch(model, val, train=False, optimizer=opt, batch_size=args.batch_size, gen=gen, device=device)
        print(f"epoch {epoch:3d} | detect F1 {m['det_F1']:.3f} (P {m['det_P']:.3f} R {m['det_R']:.3f}) "
              f"| exact {m['exact_acc']:.3f} | type|detected {m['type_acc_given_detected']:.3f}")
    print(f"[detect+type] final: detection F1 {m['det_F1']:.3f}, "
          f"exact (detect+type) {m['exact_acc']:.3f}, typing|detected {m['type_acc_given_detected']:.3f}",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
