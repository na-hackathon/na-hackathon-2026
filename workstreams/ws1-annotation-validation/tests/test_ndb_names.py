"""Item names and class labels against docs/ndb-bp-ext.dic (ndb_names.py).

Run with:  python -m pytest tests/test_ndb_names.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import gemmi
import pytest

WS1 = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WS1 / "bin"))
import ndb_names as nn  # noqa: E402

DICTIONARY = WS1 / "docs" / "ndb-bp-ext.dic"


def dictionary_items(category: str):
    """(all items, mandatory items) of a category of the dictionary."""
    prefix = f"_{category}."
    items = {}
    for item in gemmi.cif.read(str(DICTIONARY)).sole_block():
        frame = item.frame                                  # one save frame per item
        name = frame.find_value("_item.name") if frame is not None else None
        if name and gemmi.cif.as_string(name).startswith(prefix):
            items[gemmi.cif.as_string(name)[len(prefix):]] = gemmi.cif.as_string(frame.find_value("_item.mandatory_code"))
    return set(items), {i for i, m in items.items() if m == "yes"}


@pytest.mark.parametrize("category, written", [("ndb_base_pair_list", nn.LIST_ITEMS),
                                                ("ndb_base_pair_annotation", nn.ANNOTATION_ITEMS)])
def test_written_items_are_defined_and_every_mandatory_item_is_written(category, written):
    defined, mandatory = dictionary_items(category)
    assert set(written) <= defined, set(written) - defined
    assert mandatory <= set(written), mandatory - set(written)


def test_legacy_names_map_to_dictionary_names():
    assert set(nn.LEGACY_LIST_ITEMS) <= set(nn.LIST_ITEMS) and set(nn.LEGACY_ANNOTATION_ITEMS) <= set(nn.ANNOTATION_ITEMS)
    assert nn.LEGACY_LIST_ITEMS["label_comp_id_2"] == "comp_id_2" and nn.LEGACY_ANNOTATION_ITEMS["sub_class"] == "subclass"


def block_from(text: str):
    return gemmi.cif.read_string("data_x\n" + text).sole_block()


@pytest.mark.parametrize("comp, letter", [("G", "G"), ("U", "U"), ("DA", "A"), ("DT", "T"),
                                          ("PSU", "U"), ("2MG", "G"), ("OMC", "C"), ("1MA", "A"),   # gemmi table
                                          ("GTP", "G"), ("ADP", "A"), ("UMP", "U"), ("QQQ", "N")])    # nucleoside phosphates: base letter
def test_parent_base_without_a_mapping_in_the_structure(comp, letter):
    assert nn.parent_base_lookup(block_from("_cell.length_a 1\n"))(comp) == letter


def test_parent_base_from_the_structure_mapping_takes_precedence():
    block = block_from("loop_\n_pdbx_struct_mod_residue.id\n_pdbx_struct_mod_residue.label_comp_id\n"
                       "_pdbx_struct_mod_residue.parent_comp_id\n1 QQQ A\n"
                       "loop_\n_chem_comp.id\n_chem_comp.mon_nstd_parent_comp_id\nZZZ C\nPSU G\nA ?\n")
    letter = nn.parent_base_lookup(block)
    assert (letter("QQQ"), letter("ZZZ"), letter("PSU"), letter("A")) == ("A", "C", "G", "A")


@pytest.mark.parametrize("family, a, b, expected", [
    ("cWW", "G", "C", "cWW_G-C"), ("cWW", "C", "G", "cWW_G-C"),          # symmetric: bases in dictionary order
    ("cWW", "U", "A", "cWW_A-U"), ("tSS", "C", "A", "tSS_A-C"), ("cWW", "C", "U", "cWW_C-U"),
    ("tWH", "U", "A", "tWH_U-A"), ("cWS", "C", "G", "cWS_C-G"),          # asymmetric: first base has the first edge
    ("cWW", "2MG", "C", "cWW_G-C"),                                       # the example of 1EHZ (A10 2MG - A25 C)
    ("cWW", "DA", "DT", "?"), ("cWW", "G", "QQQ", "?"),                  # T and unknown parents are not in the dictionary
])
def test_pair_class(family, a, b, expected):
    assert nn.pair_class(family, a, b, nn.parent_base_lookup(block_from("_cell.length_a 1\n"))) == expected
