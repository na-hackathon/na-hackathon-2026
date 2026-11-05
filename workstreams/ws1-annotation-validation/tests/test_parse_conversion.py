"""Regression test for the base-pair exporters in bin/.

Runs the FR3D and RNApolis exporters on the sample inputs under
data/tests/{fr3d_outputs,rnapolis_outputs} (structures in data/tests/Inputs/) and compares the pairs
with those of data/tests/conversion_outputs: family, model, both residues (auth and label ids,
component, insertion code, alternate position) and symmetry operators, with their multiplicity. The
references use the item names from before docs/ndb-bp-ext.dic, so the comparison goes through
parse_basepairs, which reads both; class labels and the dictionary conformance are in
tests/test_dictionary.py.
Run with:  python -m pytest tests/test_parse_conversion.py
"""
from __future__ import annotations

import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
WS1 = Path(__file__).resolve().parents[1]
SCRIPTS = WS1 / "bin"
DATA = REPO / "data" / "tests"

INPUTS = DATA / "Inputs"
FR3D_RAW = DATA / "fr3d_outputs"
RNAPOLIS_RAW = DATA / "rnapolis_outputs"
FR3D_REF = DATA / "conversion_outputs" / "fr3d_mmcif_outputs"
RNAPOLIS_REF = DATA / "conversion_outputs" / "rnapolis_mmcif_outputs"

PDB_IDS = ["8bwt", "8vjt", "8xzn", "9e5i", "9hrf", "9sfq"]

sys.path.insert(0, str(SCRIPTS))
from parse_basepairs import parse_basepairs  # noqa: E402

FIELDS = ("family", "model") + tuple(
    f"{item}_{side}" for side in "12"
    for item in ("auth_asym_id", "auth_seq_id", "comp_id", "label_asym_id", "label_seq_id", "ins_code", "alt_id", "struct_oper_id"))


def _pairs(cif_path: Path) -> Counter:
    """The pairs of a base-pair mmCIF with their multiplicity (order-independent, ids dropped)."""
    return Counter(tuple(pair[field] for field in FIELDS) for pair in parse_basepairs(str(cif_path)))


def _run(*cmd: str) -> None:
    env_paths = f"{SCRIPTS}:{sys.path[0] if sys.path else ''}"
    subprocess.run(
        list(cmd),
        check=True,
        cwd=SCRIPTS,
        env={**__import__("os").environ, "PYTHONPATH": env_paths},
    )


@pytest.fixture(scope="module")
def workdir(tmp_path_factory) -> Path:
    return tmp_path_factory.mktemp("conversion_outputs")


@pytest.mark.parametrize("pdb_id", PDB_IDS)
def test_fr3d_conversion(pdb_id: str, workdir: Path) -> None:
    raw = FR3D_RAW / pdb_id / f"{pdb_id}_basepair.txt"
    structure = INPUTS / f"{pdb_id}.cif"
    reference = FR3D_REF / f"{pdb_id}_fr3d_basepairs.cif"
    assert raw.exists() and structure.exists() and reference.exists()

    out_dir = workdir / "fr3d_mmcif_outputs"
    out_dir.mkdir(exist_ok=True)
    json_path = out_dir / f"{pdb_id}.json"
    out_cif = out_dir / f"{pdb_id}_fr3d_basepairs.cif"

    _run(sys.executable, str(SCRIPTS / "fr3d_basepair_to_json.py"), str(raw), str(json_path))
    _run(
        sys.executable,
        str(SCRIPTS / "json_fr3d_mmcif.py"),
        str(json_path),
        str(structure),
        str(out_cif),
    )

    assert out_cif.exists(), f"parser did not produce {out_cif}"
    assert _pairs(out_cif) == _pairs(reference)


@pytest.mark.parametrize("pdb_id", PDB_IDS)
def test_rnapolis_conversion(pdb_id: str, workdir: Path) -> None:
    raw = RNAPOLIS_RAW / f"{pdb_id}.json"
    structure = INPUTS / f"{pdb_id}.cif"
    reference = RNAPOLIS_REF / f"{pdb_id}_annotated.cif"
    assert raw.exists() and structure.exists() and reference.exists()

    out_dir = workdir / "rnapolis_mmcif_outputs"
    out_dir.mkdir(exist_ok=True)
    out_cif = out_dir / f"{pdb_id}_annotated.cif"

    _run(
        sys.executable,
        str(SCRIPTS / "json_rnapolis_mmcif.py"),
        str(raw),
        str(structure),
        str(out_cif),
    )

    assert out_cif.exists(), f"parser did not produce {out_cif}"
    assert _pairs(out_cif) == _pairs(reference)
