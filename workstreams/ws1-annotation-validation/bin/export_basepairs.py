#!/usr/bin/env python3
"""Reparse a base-pair mmCIF and emit:
    <prefix>.csv   pair table (CSV)
    <prefix>.tsv   pair table (TSV, same as parse_basepairs.py stdout)
    <prefix>.bg    forgi BulgeGraph built from cWW pairs
    <prefix>.jpg   matplotlib render of that BulgeGraph

Usage: export_basepairs.py <basepairs.cif> <out_prefix>
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import gemmi

# parse_basepairs lives next door in read_write_mmcif/
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "read_write_mmcif"))
from parse_basepairs import parse_basepairs  # noqa: E402

NA_RESIDUES = {"A", "U", "G", "C", "T", "DA", "DT", "DG", "DC"}

CSV_COLUMNS = [
    "base_pair_id", "family",
    "auth_asym_id_1", "auth_seq_id_1", "comp_id_1",
    "auth_asym_id_2", "auth_seq_id_2", "comp_id_2",
]


def write_csv(pairs: list[dict], path: Path) -> None:
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CSV_COLUMNS)
        w.writeheader()
        w.writerows(pairs)


def write_tsv(pairs: list[dict], path: Path) -> None:
    with open(path, "w") as f:
        f.write("res_1\tfamily\tres_2\n")
        for p in pairs:
            f.write(
                f"{p['auth_asym_id_1']}|{p['auth_seq_id_1']}|{p['comp_id_1']}\t"
                f"{p['family']}\t"
                f"{p['auth_asym_id_2']}|{p['auth_seq_id_2']}|{p['comp_id_2']}\n"
            )


def _chain_dotbracket(chain, pairs: list[dict]) -> tuple[str, str] | None:
    """For one chain: build (sequence, dot-bracket) from its cWW pairs.
    Returns None if the chain has no nucleic-acid residues.
    """
    residues = [(r.seqid.num, r.name) for r in chain if r.name in NA_RESIDUES]
    if not residues:
        return None
    idx = {seq: i for i, (seq, _) in enumerate(residues)}
    seq_str = "".join(name[-1] for _, name in residues)  # DA -> A, U -> U
    db = ["."] * len(residues)
    for p in pairs:
        if p["family"] != "cWW":
            continue
        if p["auth_asym_id_1"] != chain.name or p["auth_asym_id_2"] != chain.name:
            continue
        i = idx.get(int(p["auth_seq_id_1"]))
        j = idx.get(int(p["auth_seq_id_2"]))
        if i is None or j is None:
            continue
        a, b = sorted((i, j))
        db[a], db[b] = "(", ")"
    return seq_str, "".join(db)


def build_bulge_graph(pairs: list[dict], cif_path: Path):
    """Build a forgi BulgeGraph from the cWW pairs of every nucleic-acid chain.
    Returns the BulgeGraph, or None when no nucleic-acid chain was found.
    """
    from forgi.graph.bulge_graph import BulgeGraph

    structure = gemmi.read_structure(str(cif_path))
    blocks = [b for b in (_chain_dotbracket(c, pairs) for c in structure[0]) if b]
    if not blocks:
        return None
    seq = "&".join(s for s, _ in blocks)
    dotbracket = "&".join(d for _, d in blocks)
    bg = BulgeGraph.from_dotbracket(dotbracket, seq)
    bg.name = cif_path.stem
    return bg


def write_forgi_bg(bg, out_path: Path) -> None:
    if bg is None:
        out_path.write_text("# no nucleic-acid chains found\n")
        return
    out_path.write_text(bg.to_bg_string())


def write_forgi_jpeg(bg, out_path: Path) -> None:
    """Render the BulgeGraph as a JPEG via forgi.visual.mplotlib."""
    if bg is None:
        # empty placeholder so downstream tools see the file
        out_path.write_bytes(b"")
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import forgi.visual.mplotlib as fvm

    fig, ax = plt.subplots(figsize=(8, 8))
    fvm.plot_rna(bg, ax=ax)
    ax.set_title(bg.name or "")
    ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, format="jpg")
    plt.close(fig)


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} <basepairs.cif> <out_prefix>")
    cif = Path(sys.argv[1])
    prefix = Path(sys.argv[2])

    pairs = parse_basepairs(str(cif))
    write_csv(pairs, prefix.with_suffix(".csv"))
    write_tsv(pairs, prefix.with_suffix(".tsv"))
    bg = build_bulge_graph(pairs, cif)
    write_forgi_bg(bg, prefix.with_suffix(".bg"))
    write_forgi_jpeg(bg, prefix.with_suffix(".jpg"))


if __name__ == "__main__":
    main()
