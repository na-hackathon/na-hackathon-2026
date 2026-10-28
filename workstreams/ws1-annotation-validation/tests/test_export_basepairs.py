"""Unit tests for the --visualize exports of bin/export_basepairs.py (one forgi graph per chain).

Needs forgi (envs/ws1.yml); skipped where it is not installed.
Run with:  python -m pytest tests/test_export_basepairs.py
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("forgi")
BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))
import export_basepairs as eb  # noqa: E402

REFS = Path(__file__).resolve().parents[3] / "data" / "tests" / "conversion_outputs"
HAIRPINS = {"A": "GGAAACC", "B": "CCUUUGG"}      # each chain pairs 1-7 and 2-6; no pair between the chains


def write_pdb(path: Path, chains: dict) -> Path:
    """Structure with one P atom per residue, numbered from 1; a chain is a list of residue names."""
    lines, serial = [], 0
    for chain, residues in chains.items():
        for number, name in enumerate(residues, start=1):
            serial += 1
            lines.append(f"ATOM  {serial:5d}  P   {name:>3s} {chain}{number:4d}"
                         f"{float(number):12.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           P")
    path.write_text("\n".join(lines) + "\nEND\n")
    return path


def cww(chain, i, j):
    return {"family": "cWW", "auth_asym_id_1": chain, "auth_seq_id_1": str(i),
            "auth_asym_id_2": chain, "auth_seq_id_2": str(j)}


PAIRS = [cww(c, i, j) for c in HAIRPINS for i, j in ((1, 7), (2, 6))]


def test_one_graph_per_chain_when_chains_are_not_paired(tmp_path):
    graphs = eb.build_bulge_graphs(PAIRS, write_pdb(tmp_path / "s.pdb", HAIRPINS))
    assert [g.name for g in graphs] == ["s_A", "s_B"]
    assert [g.to_dotbracket_string() for g in graphs] == ["((...))", "((...))"]


def test_single_chain_keeps_the_structure_name(tmp_path):
    structure = write_pdb(tmp_path / "s.pdb", {"A": HAIRPINS["A"]})
    (graph,) = eb.build_bulge_graphs([p for p in PAIRS if p["auth_asym_id_1"] == "A"], structure)
    assert graph.name == "s"


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
