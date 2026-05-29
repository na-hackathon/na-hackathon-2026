#!/usr/bin/env python3
"""
ndb_to_rnaglib.py

rnaglib cannot ingest the NDB/NAPAIR base-pair annotation mmCIF directly: its
only file path runs FR3D, which needs `_atom_site` coordinates this file does
not have. But the file already contains what an rnaglib 2.5D graph stores
(base pairs + Leontis-Westhof families), so this adapter rewrites it into an
rnaglib-native graph that rnaglib's own `load_graph` can read.

Output graph (rnaglib convention):
    nodes : "<pdbid>.<chain>.<seq>" with nt_code / chain_id / is_modified
    edges : LW-labelled, both directions
            base pairs -> e.g. cWW / tWS (reverse direction = tSW), and
            backbone   -> B53 (5'->3') and B35 (3'->5')
    graph : graph["pdbid"] set (required by rnaglib.utils.load_graph)

Usage:
    python ndb_to_rnaglib.py <input.mmcif> [-o graph.json] [--model N]
    # then, read it back through rnaglib:
    #   from rnaglib.utils import load_graph; G = load_graph("graph.json")
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import networkx as nx

from mmcif_to_pyg import parse_ndb_mmcif

STANDARD_NT = {"A", "C", "G", "U"}


def _reverse_lw(family: str) -> str:
    """Reverse an LW family for the opposite edge direction (tWS -> tSW)."""
    if len(family) == 3:  # orientation + edge1 + edge2
        return family[0] + family[2] + family[1]
    return family


def ndb_mmcif_to_rnaglib_graph(
    path: str | Path, model: int = 1, pdbid: str | None = None
) -> nx.DiGraph:
    """Build an rnaglib-compatible 2.5D graph from an NDB annotation mmCIF."""
    t = parse_ndb_mmcif(path)
    pdbid = (pdbid or t.entry_id).lower()
    bp, unp, ann = t.base_pair_list, t.unpaired_list, t.annotation
    ann_idx = {str(v): i for i, v in enumerate(ann.get("base_pair_id", []))}

    G: nx.DiGraph = nx.DiGraph()
    G.graph["pdbid"] = pdbid
    G.graph["name"] = pdbid

    def node_id(chain: str, seq) -> str:
        return f"{pdbid}.{chain}.{int(seq)}"

    def add_residue(chain, seq, comp) -> str:
        comp = str(comp)
        n = node_id(chain, seq)
        if n not in G:
            G.add_node(
                n,
                nt_code=comp if comp in STANDARD_NT else "N",
                nt=comp,
                nt_full=comp,
                chain_id=str(chain),
                index=int(seq),
                is_modified=comp not in STANDARD_NT,
            )
        return n

    # base pairs (also registers the paired residues as nodes) ---------------
    for r in range(t.n_pairs()):
        if int(bp["PDB_model_num"][r]) != model:
            continue
        n1 = add_residue(bp["label_asym_id_1"][r], bp["label_seq_id_1"][r], bp["label_comp_id_1"][r])
        n2 = add_residue(bp["label_asym_id_2"][r], bp["label_seq_id_2"][r], bp["label_comp_id_2"][r])
        a = ann_idx.get(str(bp["id"][r]))
        family = ann["l-w_family_name"][a] if a is not None else "cWW"
        G.add_edge(n1, n2, LW=family)
        G.add_edge(n2, n1, LW=_reverse_lw(family))

    # unpaired residues -> isolated nodes ------------------------------------
    for r in range(len(unp.get("id", []))):
        if int(unp["PDB_model_num"][r]) != model:
            continue
        add_residue(unp["label_asym_id"][r], unp["label_seq_id"][r], unp["label_comp_id"][r])

    # backbone edges between consecutive residues of the same chain ----------
    by_chain: dict[str, list[tuple[int, str]]] = {}
    for n, d in G.nodes(data=True):
        by_chain.setdefault(d["chain_id"], []).append((d["index"], n))
    for residues in by_chain.values():
        residues.sort()
        for (s1, a), (s2, b) in zip(residues, residues[1:]):
            if s2 == s1 + 1:
                G.add_edge(a, b, LW="B53")
                G.add_edge(b, a, LW="B35")

    return G


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=Path, help="NDB/NAPAIR base-pair mmCIF file.")
    p.add_argument("-o", "--output", type=Path,
                   help="Write an rnaglib-loadable .json graph (default: <input>.rnaglib.json).")
    p.add_argument("--model", type=int, default=1, help="PDB model number (default: 1).")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.input.is_file():
        raise SystemExit(f"error: input file not found: {args.input}")

    G = ndb_mmcif_to_rnaglib_graph(args.input, model=args.model)
    out = args.output or args.input.with_suffix(".rnaglib.json")

    # dump + read back through rnaglib to prove it is rnaglib-readable
    from rnaglib.utils import dump_json, load_graph
    dump_json(str(out), G)
    reloaded = load_graph(str(out))

    n_bp = sum(1 for *_, d in reloaded.edges(data=True) if d["LW"] not in ("B53", "B35")) // 2
    print(f"pdbid={reloaded.graph['pdbid']} model={args.model} "
          f"nodes={reloaded.number_of_nodes()} base_pairs={n_bp}", file=sys.stderr)
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
