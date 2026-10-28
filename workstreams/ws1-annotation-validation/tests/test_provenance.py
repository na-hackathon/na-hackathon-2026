"""Tests for the annotator provenance: mmcif_software.py, ws1-provenance, the exporters'
--provenance option and the submodule check of ws1-annotate.

Run with:  python -m pytest tests/test_provenance.py
"""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

import gemmi
import pytest

WS1 = Path(__file__).resolve().parents[1]
BIN = WS1 / "bin"
DATA = WS1.parents[1] / "data" / "tests"
FIXTURES = WS1 / "tests" / "fixtures"
sys.path.insert(0, str(BIN))
from mmcif_software import add_software_row  # noqa: E402

FR3D = {"name": "FR3D", "commit": "995de6cf5fb67e689ab8b6626988799ec4362f1e",
        "repository": "https://github.com/BGSU-RNA/fr3d-python.git"}


def software_rows(cif):
    table = gemmi.cif.read(str(cif)).sole_block().get_mmcif_category("_software.")
    return [{k.lower(): table[k][i] for k in table} for i in range(len(next(iter(table.values()), [])))]


def test_row_is_appended_to_an_existing_loop():
    doc = gemmi.cif.read(str(DATA / "conversion_outputs/fr3d_mmcif_outputs/8xzn_fr3d_basepairs.cif"))
    block = doc.sole_block()
    before = len(block.get_mmcif_category("_software.")["name"])
    assert add_software_row(block, FR3D, "base-pair annotation (WS1 pipeline)")
    table = gemmi.cif.read_string(doc.as_string()).sole_block().get_mmcif_category("_software.")
    assert len(table["name"]) == before + 1
    assert table["name"][:before] == ["PHENIX", "HKL-3000", "HKL-3000", "PHASER"]      # existing rows kept
    assert (table["pdbx_ordinal"][-1], table["name"][-1], table["version"][-1]) == (str(before + 1), "FR3D", FR3D["commit"])
    assert (table["type"][-1], table["location"][-1]) == ("program", FR3D["repository"])
    assert table["description"][-1] == "base-pair annotation (WS1 pipeline)"


def test_category_is_created_when_missing():
    doc = gemmi.cif.read_string("data_x\n_cell.length_a 1\n")
    assert add_software_row(doc.sole_block(), {"name": "RNApolis", "commit": None, "repository": None}, "a b")
    row = gemmi.cif.read_string(doc.as_string()).sole_block().get_mmcif_category("_software.")
    assert (row["pdbx_ordinal"], row["name"], row["version"], row["location"], row["description"]) == (
        ["1"], ["RNApolis"], [None], [None], ["a b"])


def test_category_written_as_key_value_pairs_is_left_unchanged():
    doc = gemmi.cif.read_string("data_x\n_software.name PHENIX\n_software.pdbx_ordinal 1\n")
    assert not add_software_row(doc.sole_block(), FR3D, "d")
    assert gemmi.cif.read_string(doc.as_string()).sole_block().find_value("_software.name") == "PHENIX"


@pytest.mark.parametrize("script, raw, fixture, expected", [
    ("json_fr3d_mmcif.py", None, "fr3d_provenance.json", ("FR3D", "995de6cf5fb67e689ab8b6626988799ec4362f1e")),
    ("json_rnapolis_mmcif.py", DATA / "rnapolis_outputs" / "8xzn.json", "rnapolis_provenance.json",
     ("RNApolis", "5ba519b5af6de0086d02e61640bdff56e2ef6bc5")),
])
def test_exporter_adds_the_annotator_and_keeps_the_pairs(script, raw, fixture, expected, tmp_path):
    if raw is None:                                   # FR3D: table -> json first
        raw = tmp_path / "bp.json"
        subprocess.run([sys.executable, str(BIN / "fr3d_basepair_to_json.py"),
                        str(DATA / "fr3d_outputs/8xzn/8xzn_basepair.txt"), str(raw)], check=True)
    structure = DATA / "Inputs" / "8xzn.cif"
    plain, with_prov = tmp_path / "plain.cif", tmp_path / "prov.cif"
    for out, extra in ((plain, []), (with_prov, ["--provenance", str(FIXTURES / fixture)])):
        subprocess.run([sys.executable, str(BIN / script), str(raw), str(structure), str(out), *extra],
                       check=True, capture_output=True)
    assert [r["name"] for r in software_rows(with_prov)][-1:] == [expected[0]]
    assert software_rows(with_prov)[-1]["version"] == expected[1]
    assert len(software_rows(with_prov)) == len(software_rows(plain)) + 1
    for category in ("_ndb_base_pair_list.", "_ndb_base_pair_annotation."):               # pairs unchanged
        assert (gemmi.cif.read(str(plain)).sole_block().get_mmcif_category(category)
                == gemmi.cif.read(str(with_prov)).sole_block().get_mmcif_category(category))


@pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")
def test_ws1_provenance_reads_the_submodule_checkout(tmp_path):
    sub = tmp_path / "src" / "fr3d-python"
    sub.mkdir(parents=True)
    (sub / "pyproject.toml").write_text('[project]\nname = "fr3d"\nversion = "1.2.3"\n')
    git = ["git", "-C", str(sub), "-c", "user.name=t", "-c", "user.email=t@t"]
    subprocess.run([*git[:3], "init", "-q"], check=True)
    subprocess.run([*git, "add", "pyproject.toml"], check=True)
    subprocess.run([*git, "commit", "-q", "-m", "x"], check=True)
    commit = subprocess.run([*git[:3], "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    structure = tmp_path / "s.cif"
    structure.write_text("data_x\n")
    out = tmp_path / "p.json"
    subprocess.run([sys.executable, str(BIN / "ws1-provenance"), "--tool", "fr3d", "--src", str(tmp_path / "src"),
                    "--structure", str(structure), "--command", "python annotate s.cif", "--out", str(out)], check=True)
    data = json.loads(out.read_text())
    assert data["software"]["name"] == "FR3D" and data["software"]["commit"] == commit
    assert data["software"]["package_version"] == "1.2.3" and data["software"]["repository"] is None   # no remote
    assert data["structure"] == {"name": "s.cif", "sha256": hashlib.sha256(b"data_x\n").hexdigest()}
    assert data["command"] == "python annotate s.cif" and data["tool"] == "fr3d"


def test_ws1_provenance_without_a_checkout_gives_null_values(tmp_path):
    structure = tmp_path / "s.cif"
    structure.write_text("data_x\n")
    out = tmp_path / "p.json"
    subprocess.run([sys.executable, str(BIN / "ws1-provenance"), "--tool", "rnapolis", "--src", str(tmp_path),
                    "--structure", str(structure), "--command", "c", "--out", str(out)], check=True)
    software = json.loads(out.read_text())["software"]
    assert software["name"] == "RNApolis" and software["commit"] is None and software["package_version"] is None

