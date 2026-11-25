"""Unit tests for bin/parse_basepairs.py on tiny hand-written base-pair mmCIFs (item names of
docs/ndb-bp-ext.dic and the names used before it) and on the NAPAIR example data/00001ehz.mmcif.

Run with:  python -m pytest tests/test_parse_basepairs.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
from parse_basepairs import identity_operator_id, parse_basepairs  # noqa: E402

CIF = """data_test
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.name
2 'crystal symmetry operation' 2_655
1 'identity operation' 1_555
loop_
_ndb_base_pair_list.base_pair_id
_ndb_base_pair_list.PDB_model_number
_ndb_base_pair_list.asym_id_1
_ndb_base_pair_list.seq_id_1
_ndb_base_pair_list.comp_id_1
_ndb_base_pair_list.PDB_ins_code_1
_ndb_base_pair_list.alt_id_1
_ndb_base_pair_list.struct_oper_id_1
_ndb_base_pair_list.asym_id_2
_ndb_base_pair_list.seq_id_2
_ndb_base_pair_list.comp_id_2
_ndb_base_pair_list.PDB_ins_code_2
_ndb_base_pair_list.alt_id_2
_ndb_base_pair_list.struct_oper_id_2
_ndb_base_pair_list.auth_asym_id_1
_ndb_base_pair_list.auth_seq_id_1
_ndb_base_pair_list.auth_asym_id_2
_ndb_base_pair_list.auth_seq_id_2
1 1 A 5 G ? . 1 B 20 C A B 2 P 105 Q 220
loop_
_ndb_base_pair_annotation.id
_ndb_base_pair_annotation.base_pair_id
_ndb_base_pair_annotation.l-w_family
1 1 cWW
"""

# the same pair written by an older exporter: no model, label or operator columns
OLD_CIF = """data_old
loop_
_ndb_base_pair_list.base_pair_id
_ndb_base_pair_list.comp_id_1
_ndb_base_pair_list.comp_id_2
_ndb_base_pair_list.auth_asym_id_1
_ndb_base_pair_list.auth_seq_id_1
_ndb_base_pair_list.auth_asym_id_2
_ndb_base_pair_list.auth_seq_id_2
1 G C P 105 Q 220
"""


def test_reads_full_residue_identity_with_mixed_case_tags(tmp_path):
    f = tmp_path / "bp.cif"
    f.write_text(CIF)
    (p,) = parse_basepairs(str(f))
    assert p["family"] == "cWW" and p["model"] == "1"
    assert (p["label_asym_id_1"], p["label_seq_id_1"], p["comp_id_1"]) == ("A", "5", "G")
    assert (p["label_asym_id_2"], p["label_seq_id_2"], p["comp_id_2"]) == ("B", "20", "C")
    assert (p["auth_asym_id_1"], p["auth_seq_id_1"], p["auth_asym_id_2"], p["auth_seq_id_2"]) == ("P", "105", "Q", "220")
    # '?' (no insertion code) and '.' (no alternate position) read as empty
    assert (p["ins_code_1"], p["alt_id_1"]) == ("", "")
    assert (p["ins_code_2"], p["alt_id_2"], p["struct_oper_id_1"], p["struct_oper_id_2"]) == ("A", "B", "1", "2")


def test_identity_operator_is_found_when_not_first(tmp_path):
    f = tmp_path / "bp.cif"
    f.write_text(CIF)
    assert identity_operator_id(str(f)) == "1"


def test_older_files_without_optional_columns_still_parse(tmp_path):
    f = tmp_path / "old.cif"
    f.write_text(OLD_CIF)
    (p,) = parse_basepairs(str(f))
    assert p["family"] == "?" and p["model"] == "" and p["label_seq_id_1"] == "" and p["struct_oper_id_2"] == ""
    assert identity_operator_id(str(f)) == "1"


# the same pair under the item names of docs/ndb-bp-ext.dic
DICTIONARY_CIF = """data_test
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.name
2 'crystal symmetry operation' 2_655
1 'identity operation' 1_555
loop_
_ndb_base_pair_list.id
_ndb_base_pair_list.PDB_model_num
_ndb_base_pair_list.label_asym_id_1
_ndb_base_pair_list.label_seq_id_1
_ndb_base_pair_list.label_comp_id_1
_ndb_base_pair_list.PDB_ins_code_1
_ndb_base_pair_list.label_alt_id_1
_ndb_base_pair_list.struct_oper_id_1
_ndb_base_pair_list.label_asym_id_2
_ndb_base_pair_list.label_seq_id_2
_ndb_base_pair_list.label_comp_id_2
_ndb_base_pair_list.PDB_ins_code_2
_ndb_base_pair_list.label_alt_id_2
_ndb_base_pair_list.struct_oper_id_2
_ndb_base_pair_list.auth_asym_id_1
_ndb_base_pair_list.auth_seq_id_1
_ndb_base_pair_list.auth_asym_id_2
_ndb_base_pair_list.auth_seq_id_2
1 1 A 5 G ? . 1 B 20 C A B 2 P 105 Q 220
loop_
_ndb_base_pair_annotation.id
_ndb_base_pair_annotation.base_pair_id
_ndb_base_pair_annotation.l-w_family_name
_ndb_base_pair_annotation.class
_ndb_base_pair_annotation.sub_class
1 1 cWW cWW_G-C ?
"""

NAPAIR_EXAMPLE = Path(__file__).resolve().parents[3] / "data" / "00001ehz.mmcif"


def test_dictionary_item_names_give_the_same_pair_as_the_older_names(tmp_path):
    older, current = tmp_path / "older.cif", tmp_path / "current.cif"
    older.write_text(CIF)
    current.write_text(DICTIONARY_CIF)
    assert parse_basepairs(str(current)) == parse_basepairs(str(older))


def test_reads_the_napair_example():
    pairs = parse_basepairs(str(NAPAIR_EXAMPLE))
    assert len(pairs) == 25
    (modified,) = [p for p in pairs if p["comp_id_1"] == "2MG"]          # A10 (2MG) - A25 (C), class cWW_G-C
    assert (modified["family"], modified["model"], modified["comp_id_2"]) == ("cWW", "1", "C")
    assert (modified["auth_asym_id_1"], modified["auth_seq_id_1"], modified["auth_seq_id_2"]) == ("A", "10", "25")
    assert (modified["label_seq_id_1"], modified["ins_code_1"], modified["alt_id_1"], modified["struct_oper_id_1"]) == ("10", "", "", "1")
