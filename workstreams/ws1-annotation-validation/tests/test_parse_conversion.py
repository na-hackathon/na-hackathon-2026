"""Simple unit test for the base-pairing parsers in bin/.
you can run it with pytest using python -m pytest tests/test_parse_conversion.py
Runs the FR3D and RNApolis conversion scripts on the sample inputs under
data/tests/{fr3d_outputs,rnapolis_outputs} (with structures from
data/tests/Inputs/) and compares the produced _ndb_base_pair_* rows
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import gemmi
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

BP_LIST = "_ndb_base_pair_list."
BP_ANNOT = "_ndb_base_pair_annotation."


def _rows(cif_path: Path, category: str) -> set[tuple[str, ...]]:
    """Return the rows of `category` as a set of tuples (order-independent)."""
    block = gemmi.cif.read(str(cif_path)).sole_block()
    cat = block.get_mmcif_category(category)
    if not cat:
        return set()
    key_drop = {BP_LIST: {"base_pair_id"}, BP_ANNOT: {"id", "base_pair_id"}}[category]
    keys = [k for k in cat.keys() if k not in key_drop]
    n = len(next(iter(cat.values())))
    return {tuple(str(cat[k][i]) for k in keys) for i in range(n)}


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
    assert _rows(out_cif, BP_LIST) == _rows(reference, BP_LIST)
    assert _rows(out_cif, BP_ANNOT) == _rows(reference, BP_ANNOT)


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
    assert _rows(out_cif, BP_LIST) == _rows(reference, BP_LIST)
    assert _rows(out_cif, BP_ANNOT) == _rows(reference, BP_ANNOT)
