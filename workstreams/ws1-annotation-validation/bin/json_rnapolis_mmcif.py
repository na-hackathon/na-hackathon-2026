#!/usr/bin/python3
"""Convert RNApolis annotator JSON to base-pairing mmCIF (same categories as json_fr3d_mmcif.py).

Usage: json_rnapolis_mmcif.py RNAPOLIS_OUTPUT.json structure.cif [out.cif] [--provenance provenance.json]

RNApolis annotates the first model of the asymmetric unit only: every pair is written with that
model and the identity symmetry operator, and without an alternate position. Input pairs that
are not written are saved to <out minus .cif>.skipped.tsv with the reason; a summary line is
printed to stderr.
"""

import sys, json, re
import gemmi
from collections import Counter
from pathlib import Path
from mmcif_software import add_software_row

bp_families = { "cWW": '1', "tWW": '2', "cWH": '3', "tWH": '4', "cWS": '5', "tWS": '6', "cHH": '7', "tHH": '8', "cHS": '9', "tHS": '10', "cSS": '11', "tSS": '12' }
long_names = {"W": "Watson-Crick", "H": "Hoogsteen", "S": "Sugar", "t": "trans", "c": "cis"}

processed = []

# input pairs that are not written to the output, with the reason; saved next to the output
# as <out>.skipped.tsv. Pairs listed in both directions are only counted (duplicates).
skipped = []
duplicates = 0

base_pair_id = 0
annotation_id = 0

pattern = re.compile("[ct][WHS][WHS]")

try:
  # optional --provenance <json from ws1-provenance>: the annotator is added to _software
  args = sys.argv[1:]
  provenance_file = None
  if '--provenance' in args:
    i = args.index('--provenance')
    provenance_file = args[i + 1]
    del args[i:i + 2]

  with open(args[0], 'r') as file:
    annotation = json.load(file)

  pairs = annotation["base_pairs"]

  mmcif_file = args[1]
  out_file = args[2] if len(args) > 2 else Path(mmcif_file).stem + '_rnapolis_basepairs.cif'

  # TODO read the file once: structure.make_mmcif_document() does not return the complete
  # document, so the file is read as a structure and again as a cif document
  structure = gemmi.read_structure(mmcif_file)
  doc = gemmi.cif.read(mmcif_file)
  block = doc.sole_block()

  asym_loop = block.get_mmcif_category('_struct_asym.')
  entity_of_asym = dict(zip(asym_loop.get('id', []), asym_loop.get('entity_id', [])))
  struct_oper_list_loop = block.get_mmcif_category('_pdbx_struct_oper_list.')

  # RNApolis annotates a single model within the asymmetric unit, so every pair
  # references the identity symmetry operation
  identity_id = '1'
  if (struct_oper_list_loop):
    for idx, sym_id in enumerate(struct_oper_list_loop['id']):
      # as_string removes quotes; identity is not always first (e.g. 5a39)
      if (gemmi.cif.as_string(struct_oper_list_loop['type'][idx]) == 'identity operation'):
        identity_id = str(sym_id)

  # the JSON does not give the annotated model; use the first model of the structure
  model_number = structure[0].num

  ndb_base_pair_list_loop = block.init_mmcif_loop('_ndb_base_pair_list.', ['base_pair_id', 'PDB_model_number', 'asym_id_1', 'entity_id_1', 'seq_id_1', 'comp_id_1', 'PDB_ins_code_1', 'alt_id_1', 'struct_oper_id_1', 'asym_id_2', 'entity_id_2', 'seq_id_2', 'comp_id_2', 'PDB_ins_code_2', 'alt_id_2', 'struct_oper_id_2', 'auth_asym_id_1', 'auth_seq_id_1', 'auth_asym_id_2', 'auth_seq_id_2'])
  ndb_base_pair_annotation_loop = block.init_mmcif_loop('_ndb_base_pair_annotation.', ['id', 'base_pair_id', 'orientation', 'base_1_edge', 'base_2_edge', 'l-w_family_num', 'l-w_family', 'class', 'subclass'])

except:
  print("usage:", sys.argv[0], "RNAPOLIS_OUTPUT.json", "structure.cif", "[out.cif]", "[--provenance provenance.json]", file=sys.stderr)
  raise

def get_residue_identifiers(structure, model_number, author_chain, author_resnum, author_resname, author_ins_code):
  try:
    model = next(m for m in structure if m.num == int(model_number))
  except StopIteration:
    return None, f"Model '{model_number}' not found."

  for chain in model:
    if chain.name.strip() == author_chain:
      for residue in chain:
        if (residue.seqid.num == author_resnum and residue.name == author_resname and (residue.seqid.icode if residue.seqid.icode != ' ' else ' ') == author_ins_code):
          if residue.subchain not in entity_of_asym:
            return None, f"No _struct_asym entry for subchain '{residue.subchain}'."
          return {
            'model': model.num,
            'asym_id': residue.subchain,
            'entity_id': entity_of_asym[residue.subchain],
            'seq_id': residue.label_seq,
            'comp_id': residue.name,
            'PDB_ins_code': (residue.seqid.icode if residue.seqid.icode != ' ' else '?'),
            # RNApolis does not report an alternate position per pair
            'alt_id': '.'
          }, None

  return None, f"No match found for model '{model_number}', chain '{author_chain}', residue '{author_resnum}', name '{author_resname}', insertion code '{author_ins_code}'."

def residue_key(residue):
  ic = residue.get('icode') or ''
  return f"{residue['chain']}|{residue['number']}|{ic}|{residue['name']}"

try:
  for pair in pairs:
    res1 = pair["nt1"]["auth"]
    res2 = pair["nt2"]["auth"]
    lw = pair["lw"]

    key12 = residue_key(res1) + "_" + residue_key(res2)
    if key12 in processed:
      duplicates += 1
      continue

    # same family check as the FR3D exporter: anything but the 12 [ct][WHS][WHS] classes is skipped
    if not pattern.fullmatch(lw):
      print("Skipping unsupported pairing family", lw, residue_key(res1), residue_key(res2), file=sys.stderr)
      skipped.append(('unsupported_family', residue_key(res1), lw, residue_key(res2), ''))
      continue

    # put the pair into canonical order (edge_1 <= edge_2); reverse families
    if (lw in bp_families):
      res_a, res_b = res1, res2
      fam = lw
    else:
      res_a, res_b = res2, res1
      fam = lw[0] + lw[2] + lw[1]

    insa = res_a.get('icode') if res_a.get('icode') else ' '
    insb = res_b.get('icode') if res_b.get('icode') else ' '

    resulta, error = get_residue_identifiers(structure, model_number=model_number, author_chain=res_a['chain'], author_resnum=int(res_a['number']), author_resname=res_a['name'], author_ins_code=insa)
    if error:
      print(error, file=sys.stderr)
      skipped.append(('residue_not_found', residue_key(res1), lw, residue_key(res2), error))
      continue

    resultb, error = get_residue_identifiers(structure, model_number=model_number, author_chain=res_b['chain'], author_resnum=int(res_b['number']), author_resname=res_b['name'], author_ins_code=insb)
    if error:
      print(error, file=sys.stderr)
      skipped.append(('residue_not_found', residue_key(res1), lw, residue_key(res2), error))
      continue

    # seq_id can be None for pairs to non-polymer nucleotides
    if not isinstance(resulta['seq_id'], int) or not isinstance(resultb['seq_id'], int):
      print(f"Skipping pair with invalid seq_id: {residue_key(res1)} {residue_key(res2)} (seq_id_1={resulta['seq_id']}, seq_id_2={resultb['seq_id']})", file=sys.stderr)
      skipped.append(('no_label_seq_id', residue_key(res1), lw, residue_key(res2), f"seq_id_1={resulta['seq_id']}, seq_id_2={resultb['seq_id']}"))
      continue

    # add both directions to processed
    processed.append(key12)
    processed.append(residue_key(res2) + "_" + residue_key(res1))

    base_pair_id += 1
    annotation_id += 1

    # populate _ndb_base_pair_list loop (struct_oper_id is always the identity operation)
    ndb_base_pair_list_loop.add_row([str(base_pair_id), str(model_number), resulta['asym_id'], resulta['entity_id'], str(resulta['seq_id']), resulta['comp_id'], resulta['PDB_ins_code'], resulta['alt_id'], identity_id, resultb['asym_id'], resultb['entity_id'], str(resultb['seq_id']), resultb['comp_id'], resultb['PDB_ins_code'], resultb['alt_id'], identity_id, res_a['chain'], str(res_a['number']), res_b['chain'], str(res_b['number'])])

    # populate _ndb_base_pair_annotation loop
    ndb_base_pair_annotation_loop.add_row([str(annotation_id), str(base_pair_id), long_names[fam[0]], long_names[fam[1]], long_names[fam[2]], bp_families[fam], fam, fam + '_' + res_a['name'] + '-' + res_b['name'], fam + '_' + res_a['name'] + '-' + res_b['name'] + '_1'])
except Exception as e:
  print(f"error {e=}, {type(e)=}", file=sys.stderr)
  raise

if provenance_file:
  with open(provenance_file) as f:
    if not add_software_row(block, json.load(f)["software"], "base-pair annotation (WS1 pipeline)"):
      print("[json_rnapolis_mmcif] _software is written as key-value pairs: annotator not added", file=sys.stderr)

options = gemmi.cif.WriteOptions()
options.align_pairs = 33
options.align_loops = 30

# and save the final mmcif
doc.write_file(out_file, options)

# save the pairs that were not written, with the reason
skipped_file = Path(out_file).with_suffix('.skipped.tsv')
with open(skipped_file, 'w') as f:
  f.write("reason\tresidue_1\tfamily\tresidue_2\tdetail\n")
  for row in skipped:
    f.write("\t".join(row) + "\n")

print(f"[json_rnapolis_mmcif] {len(pairs)} pairs read: {base_pair_id} written, {len(skipped)} skipped {dict(Counter(row[0] for row in skipped))}, {duplicates} duplicates of written pairs -> {skipped_file}", file=sys.stderr)
if len(pairs) != base_pair_id + len(skipped) + duplicates:
  print("[json_rnapolis_mmcif] WARNING: pairs read != written + skipped + duplicates", file=sys.stderr)

# print the document to stdout
cif_in_string = doc.as_string(options)
print(cif_in_string)
