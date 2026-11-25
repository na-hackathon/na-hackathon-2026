#!/usr/bin/env python3
"""Reparse a base-pair mmCIF and emit:
    <prefix>.csv   pair table (CSV)
    <prefix>.tsv   pair table (TSV, same as parse_basepairs.py stdout)
    <prefix>.bg    forgi BulgeGraph(s) built from cWW pairs, one per nucleic-acid chain
    <prefix>.jpg   matplotlib render of the BulgeGraph(s), one panel per chain

Usage: export_basepairs.py <basepairs.cif> <out_prefix>
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import gemmi

# parse_basepairs lives next to this script in bin/
sys.path.insert(0, str(Path(__file__).resolve().parent))
from parse_basepairs import identity_operator_id, parse_basepairs  # noqa: E402

NA_RESIDUES = {"A", "U", "G", "C", "T", "DA", "DT", "DG", "DC"}

CSV_COLUMNS = [
    "base_pair_id", "family",
    "auth_asym_id_1", "auth_seq_id_1", "comp_id_1",
    "auth_asym_id_2", "auth_seq_id_2", "comp_id_2",
]


def write_csv(pairs: list[dict], path: Path) -> None:
    with open(path, "w", newline="") as f:
        # parse_basepairs also returns label ids, model, operators...; the CSV keeps its columns
        w = csv.DictWriter(f, fieldnames=CSV_COLUMNS, extrasaction="ignore")
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


# Bracket types of the dot-bracket: nested pairs use (), pairs crossing them [] and then <>;
# forgi reads all three (a crossing pair is a pseudoknot).
BRACKETS = ("()", "[]", "<>")


def _crosses(p: tuple[int, int], q: tuple[int, int]) -> bool:
    return p[0] < q[0] < p[1] < q[1] or q[0] < p[0] < q[1] < p[1]


def layer_pairs(pairs: list[tuple[int, int]]):
    """Distribute the pairs (i, j), i < j, over the bracket types so that no two pairs of a layer cross.
    Shorter pairs are placed first: nested helices get (), pairs crossing them [] and then <>.
    Returns the layers and the pairs left out (no layer without a crossing pair)."""
    layers = [[] for _ in BRACKETS]
    left_out = []
    for pair in sorted(pairs, key=lambda p: (p[1] - p[0], p[0])):
        for layer in layers:
            if not any(_crosses(pair, other) for other in layer):
                layer.append(pair)
                break
        else:
            left_out.append(pair)
    return layers, left_out


def _nucleotide_letter(residue, nucleic_chain: bool = False) -> str | None:
    """One-letter code of a nucleotide residue; None for any other residue (ligand, ion, amino acid).
    A modified nucleotide gets the letter of its parent base (N when gemmi has none). In a nucleic-acid
    chain (nucleic_chain), a polymer residue that gemmi does not know as a nucleotide, such as a
    5'-terminal GTP, is a nucleotide too: N. Polymer residues have a label_seq, ligands do not."""
    if residue.name in NA_RESIDUES:
        return residue.name[-1]                       # DA -> A, U -> U
    info = gemmi.find_tabulated_residue(residue.name)
    if info is not None and info.is_nucleic_acid():
        return info.one_letter_code.strip().upper() or "N"
    if nucleic_chain and residue.label_seq is not None and (info is None or not (info.is_amino_acid() or info.is_water())):
        return "N"
    return None


def _chain_dotbracket(chain, pairs: list[dict]) -> tuple[str, str, int] | None:
    """For one chain: build (sequence, dot-bracket, number of cWW pairs left out) from its cWW pairs.
    Residues are the nucleotides of the chain, modified ones included, matched on auth number and
    insertion code. A residue takes part in one pair of the dot-bracket (later pairs of it are left
    out), and pairs that cannot be placed without crossing in the three bracket types are left out.
    Returns None if the chain has no nucleic-acid residues.
    """
    nucleic_chain = any(_nucleotide_letter(r) for r in chain)
    nucleotides = [(r, letter) for r in chain if (letter := _nucleotide_letter(r, nucleic_chain))]
    if not nucleotides:
        return None
    residues = [r for r, _ in nucleotides]
    idx = {(r.seqid.num, r.seqid.icode.strip()): i for i, r in enumerate(residues)}
    seq_str = "".join(letter for _, letter in nucleotides)
    placed, used, left_out = [], set(), 0
    for p in pairs:
        if p["family"] != "cWW":
            continue
        if p["auth_asym_id_1"] != chain.name or p["auth_asym_id_2"] != chain.name:
            continue
        i = idx.get((int(p["auth_seq_id_1"]), p.get("ins_code_1", "")))
        j = idx.get((int(p["auth_seq_id_2"]), p.get("ins_code_2", "")))
        if i is None or j is None or i == j:
            continue
        if i in used or j in used:
            left_out += 1
            continue
        used.update((i, j))
        placed.append((min(i, j), max(i, j)))
    layers, crossing = layer_pairs(placed)
    db = ["."] * len(residues)
    for (opening, closing), layer in zip(BRACKETS, layers):
        for i, j in layer:
            db[i], db[j] = opening, closing
    return seq_str, "".join(db), left_out + len(crossing)


def build_bulge_graphs(pairs: list[dict], cif_path: Path, identity: str = "1") -> list:
    """Build one forgi BulgeGraph per nucleic-acid chain from its cWW pairs.

    Only intrachain pairs enter a chain's dot-bracket, so the chains are never connected by
    a base pair and forgi refuses to join them with '&' ("Found two sequences not connected
    by any base-pair"). Each chain therefore gets its own graph, named <stem>_<chain>
    (just <stem> when there is a single chain). Returns [] when no nucleic-acid chain was found.

    The graphs are drawn on the first model of the structure: pairs of other models and pairs
    with a residue of a symmetry mate (operator other than `identity`, the id of the identity
    operation) are not used.
    """
    from forgi.graph.bulge_graph import BulgeGraph

    structure = gemmi.read_structure(str(cif_path))
    model = str(structure[0].num)
    pairs = [p for p in pairs
             if p.get("model", "") in ("", model)
             and all(p.get(f"struct_oper_id_{side}", "") in ("", identity) for side in ("1", "2"))]
    blocks = [(c.name, b) for c in structure[0] if (b := _chain_dotbracket(c, pairs))]
    graphs = []
    for name, (seq, dotbracket, left_out) in blocks:
        if left_out:
            print(f"[export_basepairs] chain {name}: {left_out} cWW pair(s) left out of the dot-bracket "
                  "(residue already paired, or more than three crossing layers)", file=sys.stderr)
        bg = BulgeGraph.from_dotbracket(dotbracket, seq)
        bg.name = cif_path.stem if len(blocks) == 1 else f"{cif_path.stem}_{name}"
        graphs.append(bg)
    return graphs


def write_forgi_bg(graphs: list, out_path: Path) -> None:
    if not graphs:
        out_path.write_text("# no nucleic-acid chains found\n")
        return
    out_path.write_text("\n".join(bg.to_bg_string() for bg in graphs))


def write_forgi_jpeg(graphs: list, out_path: Path) -> None:
    """Render the BulgeGraphs as one JPEG (one panel per chain) via forgi.visual.mplotlib.
    Falls back to an empty placeholder when forgi's plot_rna can't import its
    ViennaRNA dependency (env-specific; doesn't fail the pipeline).
    """
    if not graphs:
        # empty placeholder so downstream tools see the file
        out_path.write_bytes(b"")
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import forgi.visual.mplotlib as fvm

    fig, axes = plt.subplots(1, len(graphs), figsize=(8 * len(graphs), 8), squeeze=False)
    try:
        for bg, ax in zip(graphs, axes[0]):
            fvm.plot_rna(bg, ax=ax)
            ax.set_title(bg.name or "")
            ax.set_axis_off()
    except (ImportError, ModuleNotFoundError) as exc:
        print(
            f"[export_basepairs] forgi.plot_rna unavailable ({exc}); "
            "skipping JPEG render (install viennarna with python bindings to enable).",
            file=sys.stderr,
        )
        plt.close(fig)
        out_path.write_bytes(b"")
        return
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
    graphs = build_bulge_graphs(pairs, cif, identity_operator_id(str(cif)))
    write_forgi_bg(graphs, prefix.with_suffix(".bg"))
    write_forgi_jpeg(graphs, prefix.with_suffix(".jpg"))


if __name__ == "__main__":
    main()
