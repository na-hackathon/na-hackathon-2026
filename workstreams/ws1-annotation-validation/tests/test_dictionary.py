"""The base-pair categories written by the exporters validate against docs/ndb-bp-ext.dic (gemmi).

Only the _ndb_base_pair_* categories are validated: the dictionary extends mmcif_pdbx.dic and does not
know the other categories of the structure, and its links to _atom_site are not checked.
Run with:  python -m pytest tests/test_dictionary.py
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import gemmi
import pytest

WS1 = Path(__file__).resolve().parents[1]
BIN = WS1 / "bin"
DATA = WS1.parents[1] / "data"
DICTIONARY = WS1 / "docs" / "ndb-bp-ext.dic"
CATEGORIES = ("_ndb_base_pair_provenance.", "_ndb_base_pair_list.", "_ndb_base_pair_annotation.")
EXPECTED = re.compile(r"category indicated as CONTRIBUTED|parent tag of .* is absent")      # not about the values
PDB_IDS = ["8bwt", "8vjt", "8xzn", "9e5i", "9hrf", "9sfq"]

pytestmark = pytest.mark.skipif(not hasattr(gemmi.cif, "Ddl"), reason="gemmi without DDL validation")


def problems(path: Path) -> list[str]:
    """Messages of the validation of the _ndb_base_pair_* categories of a file, without the expected ones."""
    messages = []
    ddl = gemmi.cif.Ddl(lambda *args: messages.append(" ".join(map(str, args))))
    ddl.read_ddl(gemmi.cif.read(str(DICTIONARY)))
    block = gemmi.cif.read(str(path)).sole_block()
    text = []
    for category in CATEGORIES:
        table = block.find_mmcif_category(category)
        if table:
            rows = (" ".join(v if v in ("?", ".") else gemmi.cif.quote(v) for v in row) for row in table)
            text.append("loop_\n" + "\n".join(table.tags) + "\n" + "\n".join(rows))
    ddl.validate_cif(gemmi.cif.read_string(f"data_{block.name}\n" + "\n".join(text) + "\n"))
    return [m for m in messages if not EXPECTED.search(m)]


def run(script: str, *args) -> None:
    subprocess.run([sys.executable, str(BIN / script), *map(str, args)], check=True, capture_output=True)


def test_the_napair_example_is_the_baseline():
    assert problems(DATA / "00001ehz.mmcif") == []


@pytest.mark.parametrize("pdb_id", PDB_IDS)
@pytest.mark.parametrize("tool", ["fr3d", "rnapolis"])
def test_exporter_output_validates(pdb_id, tool, tmp_path):
    structure, out = DATA / "tests" / "Inputs" / f"{pdb_id}.cif", tmp_path / "out.cif"
    provenance = WS1 / "tests" / "fixtures" / f"{tool}_provenance.json"
    if tool == "fr3d":
        raw = tmp_path / "bp.json"
        run("fr3d_basepair_to_json.py", DATA / "tests" / "fr3d_outputs" / pdb_id / f"{pdb_id}_basepair.txt", raw)
        run("json_fr3d_mmcif.py", raw, structure, out, "--provenance", provenance)
    else:
        run("json_rnapolis_mmcif.py", DATA / "tests" / "rnapolis_outputs" / f"{pdb_id}.json", structure, out,
            "--provenance", provenance)
    assert problems(out) == []


def test_the_check_detects_values_and_names_outside_the_dictionary(tmp_path):
    text = (DATA / "00001ehz.mmcif").read_text()
    wrong_class = tmp_path / "wrong_class.cif"
    wrong_class.write_text(text.replace("cWW_G-C_1", "cWW_X-X_9"))
    old_name = tmp_path / "old_name.cif"
    old_name.write_text(text.replace("_ndb_base_pair_annotation.sub_class", "_ndb_base_pair_annotation.subclass"))
    assert any("cWW_X-X_9" in m for m in problems(wrong_class))
    assert any("subclass" in m for m in problems(old_name))
