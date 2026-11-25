"""Edge-case inputs for the base-pair exporters in bin/.

Each case starts from the stored 8XZN structure / RNApolis output and removes or alters one
element. Run with:  python -m pytest tests/test_exporters_robustness.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import gemmi

REPO = Path(__file__).resolve().parents[3]
BIN = Path(__file__).resolve().parents[1] / "bin"
STRUCTURE = REPO / "data" / "tests" / "Inputs" / "8xzn.cif"
RNAPOLIS_RAW = REPO / "data" / "tests" / "rnapolis_outputs" / "8xzn.json"

# a cWW pair of the stored FR3D output for 8XZN (unit ids: pdb|model|chain|comp|number[|atom|alt|icode|symmetry])
UNIT_A, UNIT_B = "8XZN|1|A|G|2", "8XZN|1|A|C|52"


def run(script, *args):
    return subprocess.run([sys.executable, str(BIN / script), *map(str, args)],
                          capture_output=True, text=True)


def without_category(tmp_path, name):
    doc = gemmi.cif.read(str(STRUCTURE))
    doc.sole_block().find_mmcif_category(name).erase()
    out = tmp_path / "structure.cif"
    doc.write_file(str(out))
    return out


def n_pairs(cif):
    return len(gemmi.cif.read(str(cif)).sole_block().find_mmcif_category("_ndb_base_pair_list."))


def fr3d_export(tmp_path, rows, structure=STRUCTURE):
    tsv, js, out = tmp_path / "bp.tsv", tmp_path / "bp.json", tmp_path / "out.cif"
    tsv.write_text("".join("\t".join(r) + "\n" for r in rows))
    assert run("fr3d_basepair_to_json.py", tsv, js).returncode == 0
    return run("json_fr3d_mmcif.py", js, structure, out), out


def test_fr3d_symmetry_code_without_oper_list_is_skipped(tmp_path):
    # used to raise NameError: struct_oper_dict only existed when the oper list did
    structure = without_category(tmp_path, "_pdbx_struct_oper_list.")
    r, out = fr3d_export(tmp_path, [(f"{UNIT_A}||||2_655", "cWW", f"{UNIT_B}||||2_655", "0")], structure)
    assert r.returncode == 0 and "NameError" not in r.stderr
    assert "not found in _pdbx_struct_oper_list" in r.stderr
    assert n_pairs(out) == 0


def test_fr3d_inconsistent_symmetry_operations_exit_with_a_message(tmp_path):
    # used to raise TypeError: sys.exit() was called with two arguments
    structure = without_category(tmp_path, "_pdbx_struct_oper_list.")
    r, _ = fr3d_export(tmp_path, [(f"{UNIT_A}||||2_655", "cWW", UNIT_B, "0")], structure)
    assert r.returncode != 0 and "TypeError" not in r.stderr
    assert "Number of FR3D symmetry operations" in r.stderr


def test_fr3d_missing_struct_asym_skips_pairs(tmp_path):
    structure = without_category(tmp_path, "_struct_asym.")
    r, out = fr3d_export(tmp_path, [(UNIT_A, "cWW", UNIT_B, "0")], structure)
    assert r.returncode == 0 and "No _struct_asym entry" in r.stderr
    assert n_pairs(out) == 0


def test_rnapolis_missing_struct_asym_skips_pairs(tmp_path):
    out = tmp_path / "out.cif"
    r = run("json_rnapolis_mmcif.py", RNAPOLIS_RAW, without_category(tmp_path, "_struct_asym."), out)
    assert r.returncode == 0 and "No _struct_asym entry" in r.stderr
    assert n_pairs(out) == 0


def test_rnapolis_skips_unsupported_family(tmp_path):
    pairs = json.loads(RNAPOLIS_RAW.read_text())["base_pairs"]
    raw, out = tmp_path / "raw.json", tmp_path / "out.cif"
    raw.write_text(json.dumps({"base_pairs": [pairs[0], {**pairs[1], "lw": "n/a"}]}))
    r = run("json_rnapolis_mmcif.py", raw, STRUCTURE, out)
    assert r.returncode == 0 and "Skipping unsupported pairing family n/a" in r.stderr
    assert n_pairs(out) == 1


def test_fr3d_to_json_reports_usage_and_bad_rows(tmp_path):
    assert "usage" in run("fr3d_basepair_to_json.py").stderr
    tsv = tmp_path / "bp.tsv"
    tsv.write_text(f"{UNIT_A}\tcWW\t{UNIT_B}\t0\n\nonly\ttwo\n")     # a blank line is fine, a short row is not
    r = run("fr3d_basepair_to_json.py", tsv, tmp_path / "bp.json")
    assert r.returncode != 0 and "bp.tsv:3:" in r.stderr
