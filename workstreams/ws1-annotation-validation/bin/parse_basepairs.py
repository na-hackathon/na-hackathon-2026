#!/usr/bin/python3
"""Read the _ndb_base_pair_list / _ndb_base_pair_annotation categories of a base-pairing mmCIF.

parse_basepairs(cif) returns one dict per pair: base_pair_id, family (Leontis-Westhof, '?' if
missing), model, and for each residue (suffix _1 / _2) auth_asym_id, auth_seq_id, comp_id,
label_asym_id, label_seq_id, ins_code, alt_id and struct_oper_id. Missing values are ''.

Usage: parse_basepairs.py structure_basepairs.cif   (prints a sorted pair table to stdout)
"""

import sys
import gemmi


def _text(value):
  """gemmi returns None for '?' and False for '.'; both mean 'no value' here."""
  return '' if value is None or value is False else str(value)


def _category(block, name):
  """Category as {lower-cased tag: column}.

  gemmi keeps the tag case written in the file (the exporters write PDB_model_number,
  PDB_ins_code_1, ...) while mmCIF tags are case-insensitive, so look tags up in lower case.
  """
  return {tag.lower(): column for tag, column in block.get_mmcif_category(name).items()}


def identity_operator_id(cif_file):
  """Id of the identity operation in _pdbx_struct_oper_list, chosen as the exporters do ('1' if absent)."""
  block = gemmi.cif.read(cif_file).sole_block()
  opers = _category(block, '_pdbx_struct_oper_list.')
  identity = '1'
  for oper_id, oper_type in zip(opers.get('id', []), opers.get('type', [])):
    if _text(oper_type) == 'identity operation':
      identity = _text(oper_id)
  return identity


def parse_basepairs(cif_file):
  block = gemmi.cif.read(cif_file).sole_block()

  bp_list = _category(block, '_ndb_base_pair_list.')
  bp_annot = _category(block, '_ndb_base_pair_annotation.')

  # map base_pair_id -> Leontis-Westhof family
  family = dict(zip(bp_annot.get('base_pair_id', []), bp_annot.get('l-w_family', [])))

  def optional(tag, i):
    # columns missing from older files read as '' instead of raising
    column = bp_list.get(tag)
    return _text(column[i]) if column else ''

  pairs = []
  ids = bp_list.get('base_pair_id', [])
  for i, pid in enumerate(ids):
    pair = {
      'base_pair_id': pid,
      'family': family.get(pid, '?'),
      'model': optional('pdb_model_number', i),
    }
    for side in ('1', '2'):
      pair[f'auth_asym_id_{side}'] = bp_list[f'auth_asym_id_{side}'][i]
      pair[f'auth_seq_id_{side}'] = bp_list[f'auth_seq_id_{side}'][i]
      pair[f'comp_id_{side}'] = bp_list[f'comp_id_{side}'][i]
      # full residue identity as written by the exporters
      pair[f'label_asym_id_{side}'] = optional(f'asym_id_{side}', i)
      pair[f'label_seq_id_{side}'] = optional(f'seq_id_{side}', i)
      pair[f'ins_code_{side}'] = optional(f'pdb_ins_code_{side}', i)
      pair[f'alt_id_{side}'] = optional(f'alt_id_{side}', i)
      pair[f'struct_oper_id_{side}'] = optional(f'struct_oper_id_{side}', i)
    pairs.append(pair)
  return pairs


if __name__ == '__main__':
  try:
    cif_file = sys.argv[1]
  except IndexError:
    sys.exit(f"usage: {sys.argv[0]} structure_basepairs.cif")

  pairs = parse_basepairs(cif_file)

  # one line per pair, sorted so order between files does not matter
  lines = sorted(
    f"{p['auth_asym_id_1']}|{p['auth_seq_id_1']}|{p['comp_id_1']}\t{p['family']}\t{p['auth_asym_id_2']}|{p['auth_seq_id_2']}|{p['comp_id_2']}"
    for p in pairs
  )
  for line in lines:
    print(line)

  print(f"# {len(pairs)} base pairs", file=sys.stderr)
