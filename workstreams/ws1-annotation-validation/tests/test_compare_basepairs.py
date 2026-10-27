"""Unit tests for bin/compare_basepairs.py (pair identity + cross-tool summary).

Pairs are built as dicts, so only the integration test at the end touches gemmi / the
stored data under data/tests. Run with:  python -m pytest tests/test_compare_basepairs.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))
import compare_basepairs as cb  # noqa: E402


def pair(a, b, fam="cWW", model="1", op1="", op2="", ic1="", ic2="", alt1="", alt2="",
         names=("G", "C")):
    """a, b = (chain, number); label ids are set equal to the auth ids for brevity."""
    return {
        "family": fam, "model": model,
        "auth_asym_id_1": a[0], "auth_seq_id_1": str(a[1]), "comp_id_1": names[0],
        "auth_asym_id_2": b[0], "auth_seq_id_2": str(b[1]), "comp_id_2": names[1],
        "label_asym_id_1": a[0], "label_seq_id_1": str(a[1]), "ins_code_1": ic1,
        "alt_id_1": alt1, "struct_oper_id_1": op1,
        "label_asym_id_2": b[0], "label_seq_id_2": str(b[1]), "ins_code_2": ic2,
        "alt_id_2": alt2, "struct_oper_id_2": op2,
    }


def keys(pairs, mode="label"):
    return {cb.pair_key(p, mode)[0] for p in pairs}


def test_key_is_order_independent():
    assert cb.pair_key(pair(("A", 1), ("A", 9)))[0] == cb.pair_key(
        pair(("A", 9), ("A", 1), names=("C", "G")))[0]


def test_symmetry_operator_separates_pairs_that_collide_on_auth_ids():
    asu = pair(("A", 3), ("A", 3))
    mate = pair(("A", 3), ("A", 3), op2="2_655")
    assert len(keys([asu, mate], "auth")) == 1      # earlier key: collision
    assert len(keys([asu, mate], "label")) == 2


def test_model_and_insertion_code_are_part_of_identity():
    assert len(keys([pair(("A", 5), ("A", 9)), pair(("A", 5), ("A", 9), model="2")])) == 2
    assert len(keys([pair(("A", 5), ("A", 9)), pair(("A", 5), ("A", 9), ic1="A")])) == 2


def test_alt_position_is_folded_and_counted_as_duplicate():
    idx = cb.index_pairs([pair(("A", 5), ("A", 9), alt1="A"), pair(("A", 5), ("A", 9), alt1="B")], "label")
    assert len(idx) == 1
    stats, _ = cb.summarize({"t": idx}, "label", compare=False)
    assert stats["duplicate_keys"] == {"t": 1} and stats["pair_counts"] == {"t": 1}


def test_auth_key_keeps_last_row_for_shared_keys():
    rows = cb.index_pairs([pair(("A", 5), ("A", 9), fam="cWW", model="1"),
                           pair(("A", 5), ("A", 9), fam="tWH", model="2")], "auth")
    assert cb.families(next(iter(rows.values())), "auth") == {"tWH"}
    rows = cb.index_pairs([pair(("A", 5), ("A", 9), fam="cWW"),
                           pair(("A", 5), ("A", 9), fam="tWH", alt1="B")], "label")
    assert cb.families(next(iter(rows.values())), "label") == {"cWW", "tWH"}


def test_choose_model():
    by_tool = {"fr3d": [pair(("A", 1), ("A", 9), model=m) for m in ("1", "2", "10")],
               "rnapolis": [pair(("A", 1), ("A", 9), model="1")]}
    assert cb.choose_model(by_tool, "first") == "1"
    assert cb.choose_model(by_tool, "all") is None
    assert cb.choose_model(by_tool, "3") == "3"
    assert cb.choose_model({**by_tool, "empty": []}, "first") == "1"   # a tool without pairs is ignored
    with pytest.raises(SystemExit):
        cb.choose_model({"a": [pair(("A", 1), ("A", 9), model="1")],
                         "b": [pair(("A", 1), ("A", 9), model="2")]}, "first")


def test_operators_and_asu_filter():
    p = pair(("A", 1), ("A", 9), op1="1", op2="2_655")
    cb.normalize_operators([p], identity="1")
    assert p["struct_oper_id_1"] == "" and p["struct_oper_id_2"] == "2_655"
    assert not cb.in_asu(p) and cb.in_asu(pair(("A", 1), ("A", 9)))


def test_summary_shared_discordant_only_and_name_mismatch():
    t1 = [pair(("A", 1), ("A", 9)),                          # shared, concordant
          pair(("A", 2), ("A", 8), fam="cWW"),               # shared, discordant
          pair(("A", 3), ("A", 7)),                          # only in t1
          pair(("A", 4), ("A", 6), names=("G", "C"))]        # shared, residue names differ
    t2 = [pair(("A", 1), ("A", 9)),
          pair(("A", 2), ("A", 8), fam="tWH"),
          pair(("A", 5), ("A", 10)),                         # only in t2
          pair(("A", 4), ("A", 6), names=("A", "U"))]
    stats, table = cb.summarize({"t1": cb.index_pairs(t1, "label"), "t2": cb.index_pairs(t2, "label")},
                                "label", compare=True)
    assert stats["total_distinct_pairs"] == 5
    assert stats["shared_by_all_tools"] == 3
    assert stats["agreed_by_all_tools"] == 2                 # (1,9) and the name-mismatch pair share a family
    assert stats["discordant_family"] == 1
    assert stats["only_in"] == {"t1": 1, "t2": 1}
    assert stats["comp_id_mismatch"] == 1
    assert len(table) == 5


def test_single_tool_has_no_comparison_fields():
    stats, table = cb.summarize({"t": cb.index_pairs([pair(("A", 1), ("A", 9))], "label")}, "label", compare=False)
    assert stats == {"pair_counts": {"t": 1}, "duplicate_keys": {"t": 0}} and table is None


# --- stored cases (data/tests/conversion_outputs) ---------------------------------------

REFS = REPO / "data" / "tests" / "conversion_outputs"
# (distinct pairs, agreed by all tools) of the earlier auth-based key, all models
LEGACY = {"8bwt": (12, 12), "8vjt": (7, 0), "8xzn": (27, 21), "9e5i": (72, 48), "9hrf": (33, 32), "9sfq": (40, 33)}


def run_compare(pdb_id, tmp_path, *flags):
    out = tmp_path / pdb_id
    subprocess.run([sys.executable, str(BIN / "compare_basepairs.py"), "--compare", "--out", str(out),
                    "--annotations",
                    str(REFS / "fr3d_mmcif_outputs" / f"{pdb_id}_fr3d_basepairs.cif"),
                    str(REFS / "rnapolis_mmcif_outputs" / f"{pdb_id}_annotated.cif"), *flags],
                   check=True, capture_output=True)
    return json.loads(out.with_suffix(".json").read_text())


@pytest.mark.parametrize("pdb_id", sorted(LEGACY))
def test_key_auth_reproduces_earlier_counts(pdb_id, tmp_path):
    report = run_compare(pdb_id, tmp_path, "--key", "auth")
    assert (report["total_distinct_pairs"], report["agreed_by_all_tools"]) == LEGACY[pdb_id]


def test_default_key_keeps_models_apart(tmp_path):
    # 8BWT: FR3D annotated 10 models (120 rows), RNApolis only the first (12 rows)
    report = run_compare("8bwt", tmp_path)
    assert report["rows"] == {"8bwt_annotated": 12, "8bwt_fr3d": 120}
    assert report["model"] == "1" and report["rows_outside_compared_model"]["8bwt_fr3d"] == 108
    assert report["total_distinct_pairs"] == 12 and report["agreed_by_all_tools"] == 12
