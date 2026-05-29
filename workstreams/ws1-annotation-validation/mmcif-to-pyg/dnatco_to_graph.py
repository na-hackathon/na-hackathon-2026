#!/usr/bin/env python3
"""
dnatco_to_graph.py

Convert DNATCO "extended" mmCIF files into rnaglib-format 2.5D graph JSON
(readable by rnaglib.utils.load_graph and by train_dnatco.py).

These files are coordinate mmCIFs that embed NDB base-pair annotations, but with
slightly different column names than the standalone NAPAIR files
(`base_pair_id`, `PDB_model_number`, `asym_id_1`, `comp_id_1`, `l-w_family`, ...).

Graph built (rnaglib convention):
    nodes : every polymer residue from _pdbx_poly_seq_scheme
            attrs nt_code (A/C/G/U or N), nt, nt_full, chain_id, index, is_modified
    edges : backbone B53 (5'->3') / B35 (3'->5')
            base pairs from _ndb_base_pair_list + _ndb_base_pair_annotation
            (LW family forward, transposed family backward: tWS<->tSW)
    graph : graph["pdbid"] set (required by rnaglib.utils.load_graph)

Usage:
    python dnatco_to_graph.py cifs/ -o graphs/        # whole directory
    python dnatco_to_graph.py cifs/1ehz.cif -o graphs/
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import gemmi
import networkx as nx

STANDARD_NT = {"A", "C", "G", "U"}


def _reverse_lw(family: str) -> str:
    return family[0] + family[2] + family[1] if len(family) == 3 else family


def _col(cat: dict, *names: str) -> list:
    """Return the first present column among `names` (tolerates schema variants)."""
    for n in names:
        if n in cat:
            return cat[n]
    return []


# representative atoms kept per nucleotide (experimental geometry)
_WANTED_ATOMS = {"C1'", "N9", "N1", "P"}


def _residue_coords(block, model) -> dict:
    """Per-residue experimental coordinates from _atom_site, keyed by (label_asym, label_seq).

    Uses the *label* numbering so keys line up with the graph nodes. Returns
    {(asym, seq): {"C1'":[x,y,z], "N9"/"N1":[...], "P":[...]}}.
    """
    cat = block.get_mmcif_category("_atom_site.")
    if not cat or "Cartn_x" not in cat:
        return {}
    asym, seq, atom = cat.get("label_asym_id", []), cat.get("label_seq_id", []), cat.get("label_atom_id", [])
    xs, ys, zs = cat["Cartn_x"], cat["Cartn_y"], cat["Cartn_z"]
    models = cat.get("pdbx_PDB_model_num", [])

    out: dict = {}
    for i in range(len(atom)):
        if model is not None and models and models[i] != model:
            continue
        if atom[i] not in _WANTED_ATOMS:
            continue
        try:
            key = (asym[i], int(seq[i]))
            out.setdefault(key, {})[atom[i]] = [float(xs[i]), float(ys[i]), float(zs[i])]
        except (TypeError, ValueError):
            continue
    return out


def dnatco_cif_to_graph(path: str | Path, pdbid: str | None = None, coords: bool = True) -> nx.DiGraph:
    block = gemmi.cif.read(str(path)).sole_block()
    pdbid = (pdbid or block.name).lower()
    ps = block.get_mmcif_category("_pdbx_poly_seq_scheme.")
    bp = block.get_mmcif_category("_ndb_base_pair_list.")
    ann = block.get_mmcif_category("_ndb_base_pair_annotation.")

    G: nx.DiGraph = nx.DiGraph()
    G.graph["pdbid"] = pdbid
    G.graph["name"] = pdbid

    def node_id(chain, seq) -> str:
        return f"{pdbid}.{chain}.{int(seq)}"

    def add_residue(chain, seq, comp) -> str | None:
        try:
            si = int(seq)
        except (TypeError, ValueError):
            return None
        comp = str(comp)
        n = node_id(chain, si)
        if n not in G:
            G.add_node(
                n,
                nt_code=comp if comp in STANDARD_NT else "N",
                nt=comp,
                nt_full=comp,
                chain_id=str(chain),
                index=si,
                is_modified=comp not in STANDARD_NT,
            )
        return n

    # --- nodes: every polymer residue ---------------------------------------
    asym = _col(ps, "asym_id")
    seqs = _col(ps, "seq_id")
    mons = _col(ps, "mon_id")
    for a, s, m in zip(asym, seqs, mons):
        add_residue(a, s, m)

    # --- base-pair edges -----------------------------------------------------
    ann_idx = {str(v): i for i, v in enumerate(_col(ann, "base_pair_id"))}
    lw_col = _col(ann, "l-w_family", "l-w_family_name")
    models = _col(bp, "PDB_model_number", "PDB_model_num")
    model = min(models, key=lambda x: int(x)) if models else None

    a1, s1c, c1 = _col(bp, "asym_id_1", "label_asym_id_1"), _col(bp, "seq_id_1", "label_seq_id_1"), _col(bp, "comp_id_1", "label_comp_id_1")
    a2, s2c, c2 = _col(bp, "asym_id_2", "label_asym_id_2"), _col(bp, "seq_id_2", "label_seq_id_2"), _col(bp, "comp_id_2", "label_comp_id_2")
    ids = _col(bp, "base_pair_id", "id")

    for r in range(len(ids)):
        if model is not None and models[r] != model:
            continue
        n1 = add_residue(a1[r], s1c[r], c1[r])  # ensure node exists even if poly_seq missed it
        n2 = add_residue(a2[r], s2c[r], c2[r])
        if n1 is None or n2 is None:
            continue
        a = ann_idx.get(str(ids[r]))
        family = lw_col[a] if (a is not None and a < len(lw_col)) else "cWW"
        family = str(family)
        G.add_edge(n1, n2, LW=family)
        G.add_edge(n2, n1, LW=_reverse_lw(family))

    # --- backbone edges between consecutive residues of a chain --------------
    by_chain: dict[str, list[tuple[int, str]]] = {}
    for n, d in G.nodes(data=True):
        by_chain.setdefault(d["chain_id"], []).append((d["index"], n))
    for residues in by_chain.values():
        residues.sort()
        for (i1, x), (i2, y) in zip(residues, residues[1:]):
            if i2 == i1 + 1:
                G.add_edge(x, y, LW="B53")
                G.add_edge(y, x, LW="B35")

    # --- attach experimental coordinates (from _atom_site) -------------------
    if coords:
        cmap = _residue_coords(block, model)
        for n, d in G.nodes(data=True):
            rec = cmap.get((d["chain_id"], d["index"]))
            if not rec:
                continue
            if "C1'" in rec:
                d["xyz_C1p"] = rec["C1'"]
            gly = rec.get("N9") or rec.get("N1")
            if gly is not None:
                d["xyz_glyN"] = gly
            if "P" in rec:
                d["xyz_P"] = rec["P"]

    return G


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("input", type=Path, help="A dnatco .cif file or a directory of them.")
    p.add_argument("-o", "--out", type=Path, default=Path("graphs"),
                   help="Output directory for .json graphs (default: graphs/).")
    p.add_argument("--no-coords", action="store_true",
                   help="Skip extracting experimental coordinates (smaller graphs).")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    from rnaglib.utils import dump_json

    args = parse_args(argv)
    files = sorted(args.input.glob("*.cif")) if args.input.is_dir() else [args.input]
    if not files:
        raise SystemExit(f"error: no .cif files at {args.input}")
    args.out.mkdir(parents=True, exist_ok=True)

    ok = fail = 0
    for f in files:
        try:
            G = dnatco_cif_to_graph(f, coords=not args.no_coords)
            if G.number_of_nodes() == 0:
                raise ValueError("empty graph")
            dump_json(str(args.out / f"{f.stem}.json"), G)
            ok += 1
        except Exception as e:  # noqa: BLE001
            fail += 1
            print(f"[graph] {f.name}: {type(e).__name__}: {e}", file=sys.stderr)
    print(f"[graph] wrote {ok} graphs ({fail} failed) -> {args.out}", file=sys.stderr)
    print(args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
