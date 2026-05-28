#!/usr/bin/env python3
"""
rna_pairs.py

Shared building blocks for the RNA base-pair models
(train_pair_interactions.py, train_contact_prediction.py, train_pair_detect_type.py).

Contents:
  - vocab / constants (NT, BACKBONE relation ids, GEOM_DIM)
  - PairData               : Data whose `pair_index` batches like `edge_index`
  - parse_graph            : one pass over an rnaglib graph -> nodes/features/coords/edges/pairs
  - geom_features          : per-pair invariant distance features
  - sample_hard_negatives  : spatially-close non-pair candidates (KD-tree)
  - PairGNN                : configurable RGCN encoder + pair MLP head
  - train_eval_classification / train_eval_binary : the shared training loops
  - confusion / per_class_prf / classification_metrics / detection_typing_metrics / binary_prf
  - predict_one / per_structure_metrics / graph_node_order : per-structure helpers
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch
import torch.nn.functional as F
from scipy.spatial import cKDTree
from torch_geometric.data import Batch, Data
from torch_geometric.nn import RGCNConv

NT = ["A", "C", "G", "U", "N"]
NT_IDX = {c: i for i, c in enumerate(NT)}
BACKBONE = {"B53": 0, "B35": 1}          # backbone relation ids (5'->3', 3'->5')
GEOM_DIM = 5                             # [dC1'/10, dglyN/10, dP/10, seqsep/100, same_chain]


class PairData(Data):
    """Data whose `pair_index` (node-index pairs) batches like `edge_index`."""

    def __inc__(self, key, value, *args, **kwargs):
        if key == "pair_index":
            return self.num_nodes
        return super().__inc__(key, value, *args, **kwargs)

    def __cat_dim__(self, key, value, *args, **kwargs):
        if key == "pair_index":
            return 1
        return super().__cat_dim__(key, value, *args, **kwargs)


def _vec(d, key):
    v = d.get(key)
    return np.array(v, dtype=float) if v is not None else None


def load_graphs(graph_dir):
    """Return (sorted file list, [networkx graph, ...]) for a directory of rnaglib JSON."""
    from rnaglib.utils import load_graph
    files = sorted(Path(graph_dir).glob("*.json"))
    if not files:
        raise SystemExit(f"error: no .json graphs in {graph_dir}")
    return files, [load_graph(str(f)) for f in files]


# --- graph parsing -----------------------------------------------------------
@dataclass
class ParsedGraph:
    nodes: list
    idx: dict
    x: torch.Tensor                       # [N, len(NT)] one-hot identity
    seq: np.ndarray                       # [N] residue numbers
    chain: list                           # [N] chain ids
    c1: list                              # [N] C1' coords or None
    gly: list                             # [N] glycosidic-N coords or None
    pp: list                              # [N] phosphate coords or None
    bb_src: list
    bb_dst: list
    bb_rel: list                          # backbone directed edges (rel in BACKBONE)
    base_pairs: list                      # [(lo_idx, hi_idx, family)], lo<hi by node index

    @property
    def num_nodes(self) -> int:
        return self.x.size(0)


def parse_graph(G) -> ParsedGraph:
    """Single pass over an rnaglib graph extracting everything the models need."""
    nodes = sorted(G.nodes, key=lambda n: (G.nodes[n]["chain_id"], G.nodes[n]["index"]))
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)

    x = torch.zeros(n, len(NT))
    seq = np.zeros(n, dtype=int)
    chain = [""] * n
    c1 = [None] * n; gly = [None] * n; pp = [None] * n
    for node in nodes:
        i, d = idx[node], G.nodes[node]
        x[i, NT_IDX.get(d["nt_code"], NT_IDX["N"])] = 1.0
        seq[i] = d["index"]; chain[i] = d["chain_id"]
        c1[i] = _vec(d, "xyz_C1p"); gly[i] = _vec(d, "xyz_glyN"); pp[i] = _vec(d, "xyz_P")

    bb_src, bb_dst, bb_rel, base_pairs = [], [], [], []
    for u, v, d in G.edges(data=True):
        lw = d["LW"]
        iu, iv = idx[u], idx[v]
        if lw in BACKBONE:
            bb_src.append(iu); bb_dst.append(iv); bb_rel.append(BACKBONE[lw])
        elif iu < iv:                     # canonical direction -> one entry per pair
            base_pairs.append((iu, iv, lw))
    return ParsedGraph(nodes, idx, x, seq, chain, c1, gly, pp, bb_src, bb_dst, bb_rel, base_pairs)


# --- geometry + negatives ----------------------------------------------------
def geom_features(a, b, seq, chain, c1, gly, pp):
    """Per-pair invariant features and the raw C1'-C1' distance (Angstrom)."""
    dc1 = float(np.linalg.norm(c1[a] - c1[b]))
    dg = float(np.linalg.norm(gly[a] - gly[b])) if gly[a] is not None and gly[b] is not None else dc1
    dp = float(np.linalg.norm(pp[a] - pp[b])) if pp[a] is not None and pp[b] is not None else dc1
    ss = abs(seq[a] - seq[b]) if chain[a] == chain[b] else 0
    return [dc1 / 10, dg / 10, dp / 10, min(ss, 100) / 100, float(chain[a] == chain[b])], dc1


def sample_hard_negatives(c1, seq, chain, positives, cutoff, n_neg, rng):
    """Spatially-close (<= cutoff A on C1'), non-pair, non-adjacent node pairs."""
    if n_neg <= 0:
        return []
    has = [i for i in range(len(c1)) if c1[i] is not None]
    tree = cKDTree(np.array([c1[i] for i in has]))
    cand = []
    for a, b in tree.query_pairs(cutoff):
        ni, nj = has[a], has[b]
        lo, hi = (ni, nj) if ni < nj else (nj, ni)
        if (lo, hi) in positives:
            continue
        if chain[lo] == chain[hi] and abs(seq[lo] - seq[hi]) <= 1:
            continue
        cand.append((lo, hi))
    if not cand:
        return []
    n = min(len(cand), n_neg)
    return [cand[k] for k in rng.choice(len(cand), size=n, replace=False)]


# --- model -------------------------------------------------------------------
class PairGNN(torch.nn.Module):
    """RGCN node encoder + MLP head over [h_i, h_j, h_i*h_j, (geom)].

    out_dim=1 -> binary logit (use BCE); out_dim=C -> class logits (use CE).
    """

    def __init__(self, out_dim: int, num_relations: int, hidden: int = 64, use_geom: bool = False):
        super().__init__()
        self.num_classes = out_dim
        self.use_geom = use_geom
        self.conv1 = RGCNConv(len(NT), hidden, num_relations)
        self.conv2 = RGCNConv(hidden, hidden, num_relations)
        in_dim = 3 * hidden + (GEOM_DIM if use_geom else 0)
        self.head = torch.nn.Sequential(
            torch.nn.Linear(in_dim, hidden), torch.nn.ReLU(),
            torch.nn.Linear(hidden, out_dim),
        )

    def forward(self, x, edge_index, edge_type, pair_index, pair_geom=None):
        h = F.relu(self.conv1(x, edge_index, edge_type))
        h = F.relu(self.conv2(h, edge_index, edge_type))
        i, j = pair_index
        parts = [h[i], h[j], h[i] * h[j]]
        if self.use_geom:
            parts.append(pair_geom)
        return self.head(torch.cat(parts, dim=-1))


def _forward(model, batch):
    return model(batch.x, batch.edge_index, batch.edge_type,
                 batch.pair_index, getattr(batch, "pair_geom", None))


# --- training loops ----------------------------------------------------------
def train_eval_classification(model, data_list, *, train, optimizer, batch_size, gen, device):
    """Cross-entropy loop; returns (confusion[C, C], mean_loss)."""
    model.train(train)
    order = (torch.randperm(len(data_list), generator=gen).tolist()
             if train else list(range(len(data_list))))
    C = model.num_classes
    conf = torch.zeros(C, C, dtype=torch.long)
    tot_loss = total = 0
    for s in range(0, len(order), batch_size):
        batch = Batch.from_data_list([data_list[i] for i in order[s:s + batch_size]]).to(device)
        with torch.set_grad_enabled(train):
            logits = _forward(model, batch)
            loss = F.cross_entropy(logits, batch.pair_y)
            if train:
                optimizer.zero_grad(); loss.backward(); optimizer.step()
        pred = logits.argmax(1).cpu(); y = batch.pair_y.cpu()
        conf += torch.bincount(y * C + pred, minlength=C * C).reshape(C, C)
        tot_loss += float(loss.detach()) * y.numel(); total += y.numel()
    return conf, tot_loss / max(total, 1)


def train_eval_binary(model, data_list, *, train, optimizer, batch_size, gen, device):
    """Binary-cross-entropy loop; returns (pred[0/1], y[0/1], mean_loss)."""
    model.train(train)
    order = (torch.randperm(len(data_list), generator=gen).tolist()
             if train else list(range(len(data_list))))
    tot_loss = total = 0
    preds, ys = [], []
    for s in range(0, len(order), batch_size):
        batch = Batch.from_data_list([data_list[i] for i in order[s:s + batch_size]]).to(device)
        with torch.set_grad_enabled(train):
            logits = _forward(model, batch).squeeze(-1)
            loss = F.binary_cross_entropy_with_logits(logits, batch.pair_y)
            if train:
                optimizer.zero_grad(); loss.backward(); optimizer.step()
        tot_loss += float(loss.detach()) * batch.pair_y.numel(); total += batch.pair_y.numel()
        preds.append((logits.detach() > 0).long().cpu()); ys.append(batch.pair_y.long().cpu())
    return torch.cat(preds), torch.cat(ys), tot_loss / max(total, 1)


@torch.no_grad()
def confusion(model, data_list, batch_size: int = 16, device: str = "cpu") -> torch.Tensor:
    """Confusion matrix over the given graphs (deterministic, no shuffle)."""
    return train_eval_classification(model, data_list, train=False, optimizer=None,
                                     batch_size=batch_size, gen=None, device=device)[0]


# --- metrics -----------------------------------------------------------------
def per_class_prf(conf: torch.Tensor):
    """Per-class precision, recall, F1, support (index aligned with class ids)."""
    conf = conf.float()
    support = conf.sum(1); predicted = conf.sum(0); correct = conf.diag()
    recall = correct / support.clamp(min=1)
    precision = correct / predicted.clamp(min=1)
    denom = precision + recall
    f1 = torch.where(denom > 0, 2 * precision * recall / denom.clamp(min=1e-12), torch.zeros_like(denom))
    return precision, recall, f1, support


def classification_metrics(conf: torch.Tensor, canon_idx: int | None = None) -> dict:
    """micro acc, balanced acc (macro recall), macro-F1, non-canonical recall."""
    precision, recall, f1, support = per_class_prf(conf)
    correct = conf.float().diag()
    present = support > 0
    micro = (correct.sum() / support.sum().clamp(min=1)).item()
    balanced = recall[present].mean().item() if present.any() else 0.0
    macro_f1 = f1[present].mean().item() if present.any() else 0.0
    nc = present.clone()
    if canon_idx is not None and canon_idx < nc.numel():
        nc[canon_idx] = False
    nc_recall = (correct[nc].sum() / support[nc].sum()).item() if nc.any() and support[nc].sum() > 0 else 0.0
    return {"micro_acc": micro, "balanced_acc": balanced,
            "macro_f1": macro_f1, "noncanonical_recall": nc_recall}


def detection_typing_metrics(conf: torch.Tensor) -> dict:
    """For the NONE(0)+families confusion: detection P/R/F1, exact, typing|detected."""
    conf = conf.float()
    true_pair = conf[1:, :].sum()
    tp = conf[1:, 1:].sum(); fn = conf[1:, 0].sum(); fp = conf[0, 1:].sum()
    detP = (tp / (tp + fp)).item() if tp + fp else 0.0
    detR = (tp / (tp + fn)).item() if tp + fn else 0.0
    detF1 = 2 * detP * detR / (detP + detR) if detP + detR else 0.0
    exact = (conf.diag()[1:].sum() / true_pair).item() if true_pair else 0.0
    typed = (conf[1:, 1:].diag().sum() / tp).item() if tp else 0.0
    return {"det_P": detP, "det_R": detR, "det_F1": detF1,
            "exact_acc": exact, "type_acc_given_detected": typed}


def binary_prf(pred: torch.Tensor, y: torch.Tensor) -> dict:
    tp = int(((pred == 1) & (y == 1)).sum()); fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum()); tn = int(((pred == 0) & (y == 0)).sum())
    prec = tp / (tp + fp) if tp + fp else 0.0
    rec = tp / (tp + fn) if tp + fn else 0.0
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0.0
    acc = (tp + tn) / max(tp + fp + fn + tn, 1)
    return {"precision": prec, "recall": rec, "f1": f1, "acc": acc}


# --- per-structure helpers ---------------------------------------------------
@torch.no_grad()
def predict_one(model, data, device: str = "cpu") -> torch.Tensor:
    """Predicted class per candidate pair (classification models), aligned with pair_index."""
    model.eval()
    return _forward(model, Batch.from_data_list([data]).to(device)).argmax(1).cpu()


@torch.no_grad()
def per_structure_metrics(model, data_list, device: str = "cpu") -> dict:
    """Per-structure detection F1 and exact (detect+type) accuracy, for box plots."""
    det_f1, exact = [], []
    for d in data_list:
        pred = predict_one(model, d, device)
        y = d.pair_y
        npos = int((y > 0).sum())
        if npos == 0:
            continue
        tp = int(((pred > 0) & (y > 0)).sum())
        fp = int(((pred > 0) & (y == 0)).sum())
        fn = int(((pred == 0) & (y > 0)).sum())
        det_f1.append(2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else 0.0)
        exact.append(float(((pred == y) & (y > 0)).sum()) / npos)
    return {"det_f1": det_f1, "exact": exact}


def graph_node_order(path):
    """Residue labels aligned to the node order used by parse_graph (chain, index)."""
    from rnaglib.utils import load_graph
    G = load_graph(str(path))
    nodes = sorted(G.nodes, key=lambda n: (G.nodes[n]["chain_id"], G.nodes[n]["index"]))
    labels = [f"{G.nodes[n]['nt_code']}{G.nodes[n]['index']}" for n in nodes]
    chains = [G.nodes[n]["chain_id"] for n in nodes]
    seqs = [G.nodes[n]["index"] for n in nodes]
    return G.graph["pdbid"], labels, chains, seqs
