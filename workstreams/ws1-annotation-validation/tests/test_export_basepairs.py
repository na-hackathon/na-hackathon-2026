"""Unit tests for the --visualize exports of bin/export_basepairs.py (one forgi graph per chain).

Needs forgi (envs/ws1.yml); skipped where it is not installed.
Run with:  python -m pytest tests/test_export_basepairs.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import gemmi
import pytest

pytest.importorskip("forgi")
BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))
import export_basepairs as eb  # noqa: E402
from parse_basepairs import identity_operator_id, parse_basepairs  # noqa: E402

REFS = Path(__file__).resolve().parents[3] / "data" / "tests" / "conversion_outputs"
HAIRPINS = {"A": "GGAAACC", "B": "CCUUUGG"}      # each chain pairs 1-7 and 2-6; no pair between the chains


def write_pdb(path: Path, chains: dict) -> Path:
    """Structure with one P atom per residue. A chain is a sequence of residue names (numbered from 1)
    or a list of (name, number, insertion code)."""
    lines, serial = [], 0
    for chain, residues in chains.items():
        for k, residue in enumerate(residues, start=1):
            name, number, icode = residue if isinstance(residue, tuple) else (residue, k, "")
            serial += 1
            lines.append(f"ATOM  {serial:5d}  P   {name:>3s} {chain}{number:4d}{icode or ' '}"
                         f"{float(k):11.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           P")
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def cww(chain, i, j, ic1="", ic2="", model="1", op1="1", op2="1"):
    """An intrachain cWW pair as returned by parse_basepairs (operator 1 = identity)."""
    return {"family": "cWW", "model": model,
            "auth_asym_id_1": chain, "auth_seq_id_1": str(i), "ins_code_1": ic1, "struct_oper_id_1": op1,
            "auth_asym_id_2": chain, "auth_seq_id_2": str(j), "ins_code_2": ic2, "struct_oper_id_2": op2}


PAIRS = [cww(c, i, j) for c in HAIRPINS for i, j in ((1, 7), (2, 6))]


def test_one_graph_per_chain_when_chains_are_not_paired(tmp_path):
    graphs = eb.build_bulge_graphs(PAIRS, write_pdb(tmp_path / "s.pdb", HAIRPINS))
    assert [g.name for g in graphs] == ["s_A", "s_B"]
    assert [g.to_dotbracket_string() for g in graphs] == ["((...))", "((...))"]


def test_single_chain_keeps_the_structure_name(tmp_path):
    structure = write_pdb(tmp_path / "s.pdb", {"A": HAIRPINS["A"]})
    (graph,) = eb.build_bulge_graphs([p for p in PAIRS if p["auth_asym_id_1"] == "A"], structure)
    assert graph.name == "s"


def test_insertion_codes_distinguish_residues_with_the_same_number(tmp_path):
    # residues 3A (stem) and 3 (loop) share the number; the pair 3A-6 closes the third stem pair
    chain = [("G", 1, ""), ("G", 2, ""), ("G", 3, "A"), ("A", 3, ""), ("A", 4, ""), ("A", 5, ""),
             ("C", 6, ""), ("C", 7, ""), ("C", 8, "")]
    pairs = [cww("A", 1, 8), cww("A", 2, 7), cww("A", 3, 6, ic1="A")]
    (graph,) = eb.build_bulge_graphs(pairs, write_pdb(tmp_path / "s.pdb", {"A": chain}))
    assert graph.to_dotbracket_string() == "(((...)))"


def test_pairs_of_other_models_and_symmetry_mates_are_not_drawn(tmp_path):
    structure = write_pdb(tmp_path / "s.pdb", {"A": HAIRPINS["A"]})
    pairs = [cww("A", 1, 7), cww("A", 2, 6),
             cww("A", 3, 5, model="2"),          # another model
             cww("A", 3, 5, op2="2"),            # a residue of a symmetry mate
             cww("A", 4, 4, op2="2")]            # a residue paired with its own symmetry mate
    (graph,) = eb.build_bulge_graphs(pairs, structure)
    assert graph.to_dotbracket_string() == "((...))"


def test_crossing_pairs_get_their_own_bracket_type(tmp_path):
    structure = write_pdb(tmp_path / "s.pdb", {"A": "GGGGCCCC"})
    (graph,) = eb.build_bulge_graphs([cww("A", 1, 5), cww("A", 3, 7)], structure)
    assert graph.to_dotbracket_string() == "(.[.).]."


def test_residue_in_two_pairs_keeps_the_first_and_reports_the_other(tmp_path, capsys):
    structure = write_pdb(tmp_path / "s.pdb", {"A": HAIRPINS["A"]})
    (graph,) = eb.build_bulge_graphs([cww("A", 1, 7), cww("A", 1, 6)], structure)
    assert graph.to_dotbracket_string() == "(.....)"
    assert "1 cWW pair(s) left out" in capsys.readouterr().err


def test_layer_pairs_leaves_out_pairs_beyond_three_layers():
    layers, left_out = eb.layer_pairs([(0, 4), (1, 5), (2, 6), (3, 7)])         # all four cross each other
    assert [len(layer) for layer in layers] == [1, 1, 1] and len(left_out) == 1
    assert eb.layer_pairs([(0, 5), (1, 4)]) == ([[(1, 4), (0, 5)], [], []], [])   # nested pairs share a layer


@pytest.mark.parametrize("name, letter", [
    ("G", "G"), ("U", "U"), ("DA", "A"), ("DT", "T"),               # standard nucleotides
    ("PSU", "U"), ("2MG", "G"), ("OMC", "C"), ("1MA", "A"),         # modified: letter of the parent base
    ("GTP", None), ("HOH", None), ("MG", None), ("ALA", None), ("QQQ", None),   # not nucleotides of the chain
])
def test_nucleotide_letter(name, letter):
    residue = gemmi.Residue()
    residue.name = name
    assert eb._nucleotide_letter(residue) == letter


def test_unknown_polymer_residue_of_a_nucleic_acid_chain_is_a_nucleotide_n():
    terminal, ligand = gemmi.Residue(), gemmi.Residue()
    terminal.name = ligand.name = "GTP"
    terminal.label_seq = 1                     # a polymer residue (5'-terminal GTP); the ligand has no label_seq
    assert eb._nucleotide_letter(terminal, nucleic_chain=True) == "N"
    assert eb._nucleotide_letter(ligand, nucleic_chain=True) is None
    assert eb._nucleotide_letter(terminal, nucleic_chain=False) is None        # not in a nucleic-acid chain
    protein = gemmi.Residue()
    protein.name, protein.label_seq = "ALA", 2
    assert eb._nucleotide_letter(protein, nucleic_chain=True) is None


@pytest.mark.parametrize("pdb_id", ["8bwt", "8vjt", "8xzn", "9e5i", "9hrf", "9sfq"])
def test_every_cww_pair_of_the_stored_rnapolis_files_is_drawn(pdb_id):
    cif = REFS / "rnapolis_mmcif_outputs" / f"{pdb_id}_annotated.cif"
    pairs = parse_basepairs(str(cif))
    graphs = eb.build_bulge_graphs(pairs, cif, identity_operator_id(str(cif)))
    drawn = sum(c in ")]>" for g in graphs for c in g.to_dotbracket_string())
    assert drawn == sum(p["family"] == "cWW" for p in pairs)


def test_modified_nucleotides_are_drawn_and_ligands_are_not(tmp_path):
    # PSU (parent U) is the third stem pair; a GTP ligand follows the chain and must not enter it
    chain = ["G", "G", "PSU", "A", "A", "A", "A", "C", "C", "GTP"]
    pairs = [cww("A", 1, 9), cww("A", 2, 8), cww("A", 3, 7)]
    (graph,) = eb.build_bulge_graphs(pairs, write_pdb(tmp_path / "s.pdb", {"A": chain}))
    assert graph.to_dotbracket_string() == "(((...)))"
    assert str(graph.seq) == "GGUAAAACC"


def test_no_nucleic_acid_chain_gives_no_graph(tmp_path):
    structure = write_pdb(tmp_path / "s.pdb", {"A": ["ALA", "GLY", "SER"]})
    assert eb.build_bulge_graphs([], structure) == []


def test_bg_file_holds_every_graph(tmp_path):
    out = tmp_path / "s.bg"
    eb.write_forgi_bg(eb.build_bulge_graphs(PAIRS, write_pdb(tmp_path / "s.pdb", HAIRPINS)), out)
    assert [l for l in out.read_text().splitlines() if l.startswith("name ")] == ["name s_A", "name s_B"]
    eb.write_forgi_bg([], out)
    assert out.read_text().startswith("# no nucleic-acid chains found")


def test_jpeg_is_written_for_several_chains(tmp_path):
    pytest.importorskip("RNA")      # forgi.plot_rna needs the ViennaRNA python module
    out = tmp_path / "s.jpg"
    eb.write_forgi_jpeg(eb.build_bulge_graphs(PAIRS, write_pdb(tmp_path / "s.pdb", HAIRPINS)), out)
    assert out.stat().st_size > 0


# stored base-pair files whose chains are not paired to each other (8VJT, 9E5I) or a single chain (8XZN)
@pytest.mark.parametrize("pdb_id, n_graphs", [("8vjt", 2), ("9e5i", 2), ("8xzn", 1)])
def test_cli_on_stored_cases(pdb_id, n_graphs, tmp_path):
    cif = REFS / "fr3d_mmcif_outputs" / f"{pdb_id}_fr3d_basepairs.cif"
    subprocess.run([sys.executable, str(BIN / "export_basepairs.py"), str(cif), str(tmp_path / pdb_id)],
                   check=True, capture_output=True)
    assert (tmp_path / f"{pdb_id}.bg").read_text().count("\nname ") + 1 == n_graphs
    assert (tmp_path / f"{pdb_id}.csv").exists() and (tmp_path / f"{pdb_id}.tsv").exists()
