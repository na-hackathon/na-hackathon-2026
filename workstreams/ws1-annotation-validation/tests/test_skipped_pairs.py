"""The exporters save the input pairs they did not write (<out>.skipped.tsv) and account for every row.

Run with:  python -m pytest tests/test_skipped_pairs.py
"""
from __future__ import annotations

import csv
import json
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
BIN = Path(__file__).resolve().parents[1] / "bin"
DATA = REPO / "data" / "tests"


def run(script, *args):
    r = subprocess.run([sys.executable, str(BIN / script), *map(str, args)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stderr


def skipped_rows(out_cif):
    with open(Path(out_cif).with_suffix(".skipped.tsv"), newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fr3d(tmp_path, rows, structure):
    tsv, js, out = tmp_path / "bp.tsv", tmp_path / "bp.json", tmp_path / "fr3d_basepairs.cif"
    tsv.write_text("".join("\t".join(r) + "\n" for r in rows))
    run("fr3d_basepair_to_json.py", tsv, js)
    return run("json_fr3d_mmcif.py", js, structure, out), out


def test_fr3d_every_row_is_written_skipped_or_a_duplicate(tmp_path):
    a, b, c, d = "8XZN|1|A|G|2", "8XZN|1|A|C|52", "8XZN|1|A|G|3", "8XZN|1|A|C|51"
    stderr, out = fr3d(tmp_path, [
        (a, "cWW", b, "0"),                       # written
        (b, "cWW", a, "0"),                       # same pair, other direction: duplicate
        (c, "ncWW", d, "0"),                      # near pair: unsupported family
        ("8XZN|1|Z|G|2", "cWW", b, "0"),          # residue that is not in the structure
    ], DATA / "Inputs" / "8xzn.cif")
    assert "4 rows read: 1 pairs written, 2 skipped" in stderr and "1 duplicates" in stderr
    assert "WARNING" not in stderr
    rows = skipped_rows(out)
    assert [r["reason"] for r in rows] == ["unsupported_family", "residue_not_found"]
    assert rows[0]["unit_id_1"] == c and rows[0]["family"] == "ncWW" and rows[0]["unit_id_2"] == d
    assert "chain 'Z'" in rows[1]["detail"]


def test_rnapolis_every_pair_is_written_skipped_or_a_duplicate(tmp_path):
    pairs = json.loads((DATA / "rnapolis_outputs" / "8xzn.json").read_text())["base_pairs"]
    missing = json.loads(json.dumps(pairs[2]))
    missing["nt1"]["auth"]["chain"] = "Z"
    reverse = {**pairs[0], "nt1": pairs[0]["nt2"], "nt2": pairs[0]["nt1"]}
    raw, out = tmp_path / "raw.json", tmp_path / "rnapolis_basepairs.cif"
    raw.write_text(json.dumps({"base_pairs": [pairs[0], reverse, {**pairs[1], "lw": "n/a"}, missing]}))
    stderr = run("json_rnapolis_mmcif.py", raw, DATA / "Inputs" / "8xzn.cif", out)
    assert "4 pairs read: 1 written, 2 skipped" in stderr and "1 duplicates" in stderr
    assert "WARNING" not in stderr
    assert [r["reason"] for r in skipped_rows(out)] == ["unsupported_family", "residue_not_found"]


def test_nothing_skipped_gives_a_header_only_file(tmp_path):
    raw = DATA / "fr3d_outputs" / "8xzn" / "8xzn_basepair.txt"
    js, out = tmp_path / "bp.json", tmp_path / "fr3d_basepairs.cif"
    run("fr3d_basepair_to_json.py", raw, js)
    stderr = run("json_fr3d_mmcif.py", js, DATA / "Inputs" / "8xzn.cif", out)
    assert "46 rows read: 23 pairs written, 0 skipped" in stderr
    assert skipped_rows(out) == []


# exclusions on the stored cases (they set the denominator of any cross-tool comparison)
@pytest.mark.parametrize("pdb_id, tool, expected", [
    ("8vjt", "fr3d", {"alt_position_mismatch": 36}),     # 36 of 70 rows: pairs between alternate conformers A/B
    ("9hrf", "rnapolis", {"no_label_seq_id": 2}),        # pairs to nucleotides without a label_seq_id
])
def test_stored_case_exclusions(pdb_id, tool, expected, tmp_path):
    structure = DATA / "Inputs" / f"{pdb_id}.cif"
    out = tmp_path / f"{tool}_basepairs.cif"
    if tool == "fr3d":
        js = tmp_path / "bp.json"
        run("fr3d_basepair_to_json.py", DATA / "fr3d_outputs" / pdb_id / f"{pdb_id}_basepair.txt", js)
        run("json_fr3d_mmcif.py", js, structure, out)
    else:
        run("json_rnapolis_mmcif.py", DATA / "rnapolis_outputs" / f"{pdb_id}.json", structure, out)
    assert Counter(r["reason"] for r in skipped_rows(out)) == expected
