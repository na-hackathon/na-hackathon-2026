#!/usr/bin/python3

import sys
import gemmi


def parse_basepairs(cif_file):
  block = gemmi.cif.read(cif_file).sole_block()

  bp_list = block.get_mmcif_category('_ndb_base_pair_list.')
  bp_annot = block.get_mmcif_category('_ndb_base_pair_annotation.')

  # map base_pair_id -> Leontis-Westhof family
  family = dict(zip(bp_annot.get('base_pair_id', []), bp_annot.get('l-w_family', [])))

  pairs = []
  ids = bp_list.get('base_pair_id', [])
  for i, pid in enumerate(ids):
    pairs.append({
      'base_pair_id': pid,
      'family': family.get(pid, '?'),
      'auth_asym_id_1': bp_list['auth_asym_id_1'][i],
      'auth_seq_id_1': bp_list['auth_seq_id_1'][i],
      'comp_id_1': bp_list['comp_id_1'][i],
      'auth_asym_id_2': bp_list['auth_asym_id_2'][i],
      'auth_seq_id_2': bp_list['auth_seq_id_2'][i],
      'comp_id_2': bp_list['comp_id_2'][i],
    })
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
