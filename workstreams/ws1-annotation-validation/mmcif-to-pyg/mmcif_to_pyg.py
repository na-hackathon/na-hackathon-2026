#!/usr/bin/env python3
"""
mmcif_to_pyg.py

Convert an NDB / NAPAIR base-pair-annotation mmCIF file into PyTorch Geometric
graphs.

Each nucleotide residue becomes a node; each annotated base pair becomes an
(undirected) edge. The input files carry several PDB models (conformations of
the same molecule), so one graph is produced per model.

Input categories used (see https://ndb.rutgers.edu):
    _ndb_base_pair_list        base pairs -> edges (per model)
    _ndb_base_pair_annotation  Leontis-Westhof family etc. (keyed by base_pair_id)
    _ndb_base_pair_validation  napair_rmsd / napasco metrics (keyed by base_pair_id)
    _ndb_base_unpaired_list    unpaired residues -> isolated nodes (per model)

Node features (Data.x): one-hot of the residue comp_id (A/C/G/U + modified bases).
Edge features (Data.edge_attr), per directed edge:
    [0] orientation == cis
    [1] orientation == trans
    [2] napair_rmsd            (0.0 if missing)
    [3] napasco_metric         (0.0 if missing)
    [4] napasco_annotation == preferred
    [5] napasco_annotation == of_concern
Leontis-Westhof family is stored separately as integer ids in Data.edge_lw_family.

Usage:
    python mmcif_to_pyg.py <input.mmcif> [-o graphs.pt] [--model N]
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
import torch
from torch_geometric.data import Data

# Order of the continuous/one-hot edge feature columns in Data.edge_attr.
EDGE_ATTR_COLUMNS = (
    "orient_cis",
    "orient_trans",
    "napair_rmsd",
    "napasco_metric",
    "ann_preferred",
    "ann_of_concern",
)


def log(msg: str) -> None:
    print(f"[mmcif_to_pyg] {msg}", file=sys.stderr)


# A residue is uniquely identified within a model by chain + seq id (+ ins/alt).
ResKey = tuple[str, int, str, str]


def _res_key(asym_id, seq_id, ins_code, alt_id) -> ResKey:
    # gemmi maps the mmCIF '.' / '?' placeholders to False / None.
    ins = "" if ins_code in (None, False) else str(ins_code)
    alt = "" if alt_id in (None, False) else str(alt_id)
    return (str(asym_id), int(seq_id), ins, alt)


def _to_float(value) -> float:
    if value in (None, False, "", "?", "."):
        return 0.0
    return float(value)


@dataclass
class NDBTables:
    """Raw column-oriented view of the NDB categories in one mmCIF block."""

    entry_id: str
    base_pair_list: dict[str, list]
    annotation: dict[str, list]
    validation: dict[str, list]
    unpaired_list: dict[str, list]

    def n_pairs(self) -> int:
        return len(self.base_pair_list.get("id", []))


def parse_ndb_mmcif(path: str | Path) -> NDBTables:
    """Read the NDB base-pair categories from an mmCIF file with gemmi."""
    doc = gemmi.cif.read(str(path))
    block = doc.sole_block()

    def cat(name: str) -> dict[str, list]:
        # Trailing '.' selects the whole category; returns {tag: [values...]}.
        return block.get_mmcif_category(name + ".")

    return NDBTables(
        entry_id=block.name,
        base_pair_list=cat("_ndb_base_pair_list"),
        annotation=cat("_ndb_base_pair_annotation"),
        validation=cat("_ndb_base_pair_validation"),
        unpaired_list=cat("_ndb_base_unpaired_list"),
    )


@dataclass
class Vocab:
    """Stable string<->index mapping shared across every model in a file."""

    tokens: list[str] = field(default_factory=list)
    index: dict[str, int] = field(default_factory=dict)

    def fit(self, values) -> "Vocab":
        for v in values:
            if v in (None, False):
                continue
            v = str(v)
            if v not in self.index:
                self.index[v] = len(self.tokens)
                self.tokens.append(v)
        return self

    def get(self, value) -> int:
        return self.index.get(str(value), -1)

    def __len__(self) -> int:
        return len(self.tokens)


def _build_vocabs(t: NDBTables) -> tuple[Vocab, Vocab]:
    """Vocabulary of residue types and of Leontis-Westhof families."""
    comp_values = list(t.base_pair_list.get("label_comp_id_1", []))
    comp_values += list(t.base_pair_list.get("label_comp_id_2", []))
    comp_values += list(t.unpaired_list.get("label_comp_id", []))
    nt_vocab = Vocab().fit(sorted(set(str(v) for v in comp_values if v not in (None, False))))

    lw_vocab = Vocab().fit(
        sorted(set(str(v) for v in t.annotation.get("l-w_family_name", []) if v not in (None, False)))
    )
    return nt_vocab, lw_vocab


def build_graphs(
    t: NDBTables,
    nt_vocab: Vocab | None = None,
    lw_vocab: Vocab | None = None,
) -> list[Data]:
    """Build one PyTorch Geometric Data graph per PDB model in the file."""
    if nt_vocab is None or lw_vocab is None:
        nt_vocab, lw_vocab = _build_vocabs(t)

    # Annotation / validation rows are keyed by the base_pair_list `id`.
    ann_by_pair = _index_by(t.annotation, "base_pair_id")
    val_by_pair = _index_by(t.validation, "base_pair_id")

    bp = t.base_pair_list
    unp = t.unpaired_list

    # Collect every (model, residue) and every (model, pair) up front.
    model_nums = sorted(
        {int(m) for m in bp.get("PDB_model_num", [])}
        | {int(m) for m in unp.get("PDB_model_num", [])}
    )

    graphs: list[Data] = []
    for model in model_nums:
        graphs.append(
            _build_one_model(t, model, nt_vocab, lw_vocab, ann_by_pair, val_by_pair)
        )
    return graphs


def _index_by(table: dict[str, list], key: str) -> dict[str, int]:
    """Map a key column value -> row index for fast per-pair lookup."""
    return {str(v): i for i, v in enumerate(table.get(key, []))}


def _build_one_model(
    t: NDBTables,
    model: int,
    nt_vocab: Vocab,
    lw_vocab: Vocab,
    ann_by_pair: dict[str, int],
    val_by_pair: dict[str, int],
) -> Data:
    bp = t.base_pair_list
    unp = t.unpaired_list

    # --- gather nodes for this model ---------------------------------------
    nodes: dict[ResKey, dict] = {}

    def add_node(asym, seq, comp, ins, alt, paired: bool) -> ResKey:
        key = _res_key(asym, seq, ins, alt)
        node = nodes.get(key)
        if node is None:
            nodes[key] = {
                "key": key,
                "asym_id": key[0],
                "seq_id": key[1],
                "comp_id": str(comp),
                "is_paired": paired,
            }
        elif paired:
            node["is_paired"] = True
        return key

    n_bp = t.n_pairs()
    pair_rows: list[tuple[ResKey, ResKey, str]] = []  # (key1, key2, base_pair_id)
    for r in range(n_bp):
        if int(bp["PDB_model_num"][r]) != model:
            continue
        k1 = add_node(
            bp["label_asym_id_1"][r], bp["label_seq_id_1"][r], bp["label_comp_id_1"][r],
            bp["PDB_ins_code_1"][r], bp["label_alt_id_1"][r], paired=True,
        )
        k2 = add_node(
            bp["label_asym_id_2"][r], bp["label_seq_id_2"][r], bp["label_comp_id_2"][r],
            bp["PDB_ins_code_2"][r], bp["label_alt_id_2"][r], paired=True,
        )
        pair_rows.append((k1, k2, str(bp["id"][r])))

    n_unp = len(unp.get("id", []))
    for r in range(n_unp):
        if int(unp["PDB_model_num"][r]) != model:
            continue
        add_node(
            unp["label_asym_id"][r], unp["label_seq_id"][r], unp["label_comp_id"][r],
            unp["PDB_ins_code"][r], unp["label_alt_id"][r], paired=False,
        )

    # Canonical node ordering: by chain then sequence position.
    ordered = sorted(nodes.values(), key=lambda n: (n["asym_id"], n["seq_id"]))
    full_idx = {n["key"]: i for i, n in enumerate(ordered)}

    num_nodes = len(ordered)
    x = torch.zeros((num_nodes, max(len(nt_vocab), 1)), dtype=torch.float)
    node_comp_id, node_asym_id = [], []
    node_seq_id = torch.empty(num_nodes, dtype=torch.long)
    node_is_paired = torch.zeros(num_nodes, dtype=torch.bool)
    for i, n in enumerate(ordered):
        ci = nt_vocab.get(n["comp_id"])
        if ci >= 0:
            x[i, ci] = 1.0
        node_comp_id.append(n["comp_id"])
        node_asym_id.append(n["asym_id"])
        node_seq_id[i] = n["seq_id"]
        node_is_paired[i] = n["is_paired"]

    # --- gather edges for this model ---------------------------------------
    src, dst = [], []
    edge_attr_rows: list[list[float]] = []
    edge_lw: list[int] = []
    edge_orientation: list[str] = []
    edge_base_edges: list[tuple[str, str]] = []
    edge_class: list[str] = []
    edge_annotation: list[str] = []

    for k1, k2, pair_id in pair_rows:
        i, j = full_idx[k1], full_idx[k2]

        a = ann_by_pair.get(pair_id)
        orientation = t.annotation["orientation"][a] if a is not None else None
        lw_name = t.annotation["l-w_family_name"][a] if a is not None else None
        b1_edge = t.annotation["base_1_edge"][a] if a is not None else None
        b2_edge = t.annotation["base_2_edge"][a] if a is not None else None
        cls = t.annotation["class"][a] if a is not None else None

        v = val_by_pair.get(pair_id)
        rmsd = _to_float(t.validation["napair_rmsd"][v]) if v is not None else 0.0
        napasco = _to_float(t.validation["napasco_metric"][v]) if v is not None else 0.0
        annot = t.validation["napasco_annotation"][v] if v is not None else None

        feat = [
            1.0 if orientation == "cis" else 0.0,
            1.0 if orientation == "trans" else 0.0,
            rmsd,
            napasco,
            1.0 if annot == "preferred" else 0.0,
            1.0 if annot == "of_concern" else 0.0,
        ]
        lw_id = lw_vocab.get(lw_name) if lw_name is not None else -1

        # Undirected: store both directions with identical attributes.
        for a_idx, b_idx in ((i, j), (j, i)):
            src.append(a_idx)
            dst.append(b_idx)
            edge_attr_rows.append(feat)
            edge_lw.append(lw_id)
            edge_orientation.append(str(orientation))
            edge_base_edges.append((str(b1_edge), str(b2_edge)))
            edge_class.append(str(cls))
            edge_annotation.append(str(annot))

    edge_index = (
        torch.tensor([src, dst], dtype=torch.long)
        if src else torch.empty((2, 0), dtype=torch.long)
    )
    edge_attr = (
        torch.tensor(edge_attr_rows, dtype=torch.float)
        if edge_attr_rows else torch.empty((0, len(EDGE_ATTR_COLUMNS)), dtype=torch.float)
    )
    edge_lw_family = torch.tensor(edge_lw, dtype=torch.long)

    data = Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
    data.num_nodes = num_nodes
    data.edge_lw_family = edge_lw_family
    # Plain-python metadata (not collated as node/edge tensors).
    data.entry_id = t.entry_id
    data.model_num = model
    data.node_comp_id = node_comp_id
    data.node_asym_id = node_asym_id
    data.node_seq_id = node_seq_id
    data.node_is_paired = node_is_paired
    data.edge_orientation = edge_orientation
    data.edge_base_edges = edge_base_edges
    data.edge_class = edge_class
    data.edge_annotation = edge_annotation
    data.nt_vocab = nt_vocab.tokens
    data.lw_vocab = lw_vocab.tokens
    data.edge_attr_columns = list(EDGE_ATTR_COLUMNS)
    return data


def sequence(data: Data) -> str:
    """Residue sequence (5'->3') in node order. Modified bases are bracketed."""
    parts = [c if len(c) == 1 else f"[{c}]" for c in data.node_comp_id]
    return "".join(parts)


def mmcif_to_pyg(path: str | Path) -> list[Data]:
    """End-to-end: parse an NDB mmCIF file into a list of per-model graphs."""
    tables = parse_ndb_mmcif(path)
    return build_graphs(tables)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=Path, help="NDB/NAPAIR base-pair mmCIF file.")
    p.add_argument("-o", "--output", type=Path,
                   help="Write the list of graphs to this .pt file (torch.save).")
    p.add_argument("--model", type=int,
                   help="Print only the graph for this PDB model number.")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.input.is_file():
        raise SystemExit(f"error: input file not found: {args.input}")

    graphs = mmcif_to_pyg(args.input)
    log(f"{graphs[0].entry_id}: {len(graphs)} model(s)")

    shown = graphs if args.model is None else [g for g in graphs if g.model_num == args.model]
    for g in shown:
        n_edges = g.edge_index.size(1) // 2  # undirected pairs
        log(f"  model {g.model_num}: {g.num_nodes} nodes, {n_edges} base pairs")
        log(f"    sequence: {sequence(g)}")

    if args.output is not None:
        torch.save(graphs, args.output)
        log(f"saved {len(graphs)} graph(s) -> {args.output}")
        print(args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
