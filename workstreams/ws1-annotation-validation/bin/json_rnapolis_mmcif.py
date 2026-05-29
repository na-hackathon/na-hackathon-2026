#!/usr/bin/python3

# Convert RNApolis annotator JSON output to the same base-pairing mmCIF


import sys, json, re
from purine_pyrimidine import *
import gemmi
from pathlib import Path

bp_families = { "cWW": '1', "tWW": '2', "cWH": '3', "tWH": '4', "cWS": '5', "tWS": '6', "cHH": '7', "tHH": '8', "cHS": '9', "tHS": '10', "cSS": '11', "tSS": '12' }
long_names = {"W": "Watson-Crick", "H": "Hoogsteen", "S": "Sugar", "t": "trans", "c": "cis"}

complementary = [ ("A", "U"), ("G", "C"), ("DA", "DT"), ("DG", "DC") ]

processed = []

base_pair_id = 0
annotation_id = 0

pattern = re.compile("[ct][WHS][WHS]")

try:
  with open(sys.argv[1], 'r') as file:
    annotation = json.load(file)

  pairs = annotation["base_pairs"]

  mmcif_file = sys.argv[2]
  out_file = sys.argv[3] if len(sys.argv) > 3 else Path(mmcif_file).stem + '_rnapolis_basepairs.cif'

  # TODO replace structure with something else so we can drop the double reading
  structure = gemmi.read_structure(mmcif_file)

  # there is some error getting correct document from structure, so we read it twice
  # BUG doc = structure.make_mmcif_document()
  doc = gemmi.cif.read(mmcif_file)
  block = doc.sole_block()

  asym_loop = block.get_mmcif_category('_struct_asym.')
  struct_oper_list_loop = block.get_mmcif_category('_pdbx_struct_oper_list.')

  # rnapolis annotates a single model within the asymmetric unit, so every pair
  # references the identity symmetry operation
  identity_id = '1'
  if (struct_oper_list_loop):
    for idx, sym_id in enumerate(struct_oper_list_loop['id']):
      # as_string removes quotes; identity is not always first (e.g. 5a39)
      if (gemmi.cif.as_string(struct_oper_list_loop['type'][idx]) == 'identity operation'):
        identity_id = str(sym_id)

  # the model number rnapolis annotated is not in the JSON; use the structure's
  # first model (commonly 1)
  model_number = structure[0].num

  ndb_base_pair_list_loop = block.init_mmcif_loop('_ndb_base_pair_list.', ['base_pair_id', 'PDB_model_number', 'asym_id_1', 'entity_id_1', 'seq_id_1', 'comp_id_1', 'PDB_ins_code_1', 'alt_id_1', 'struct_oper_id_1', 'asym_id_2', 'entity_id_2', 'seq_id_2', 'comp_id_2', 'PDB_ins_code_2', 'alt_id_2', 'struct_oper_id_2', 'auth_asym_id_1', 'auth_seq_id_1', 'auth_asym_id_2', 'auth_seq_id_2'])
  ndb_base_pair_annotation_loop = block.init_mmcif_loop('_ndb_base_pair_annotation.', ['id', 'base_pair_id', 'orientation', 'base_1_edge', 'base_2_edge', 'l-w_family_num', 'l-w_family', 'class', 'subclass'])

except:
  print("usage:", sys.argv[0], "RNAPOLIS_OUTPUT.json", "structure.cif", "[out.cif]", file=sys.stderr)
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
          return {
            'model': model.num,
            'asym_id': residue.subchain,
            'entity_id': asym_loop["entity_id"][asym_loop["id"].index(residue.subchain)],
            'seq_id': residue.label_seq,
            'comp_id': residue.name,
            'PDB_ins_code': (residue.seqid.icode if residue.seqid.icode != ' ' else '?'),
            # rnapolis does not report altloc per pair
            'alt_id': '.'
          }, None

  return None, f"No match found for model '{model_number}', chain '{author_chain}', residue '{author_resnum}', name '{author_resname}', insertion code '{author_ins_code}'."

def auth_key(auth):
  ic = auth.get('icode') or ''
  return f"{auth['chain']}|{auth['number']}|{ic}|{auth['name']}"

try:
  for pair in pairs:
    auth1 = pair["nt1"]["auth"]
    auth2 = pair["nt2"]["auth"]
    lw = pair["lw"]

    key12 = auth_key(auth1) + "_" + auth_key(auth2)
    if key12 in processed:
      continue


    # put the pair into canonical order (edge_1 <= edge_2); reverse families
    if (lw in bp_families):
      autha, authb = auth1, auth2
      fam = lw
    else:
      autha, authb = auth2, auth1
      fam = lw[0] + lw[2] + lw[1]

    insa = autha.get('icode') if autha.get('icode') else ' '
    insb = authb.get('icode') if authb.get('icode') else ' '

    resulta, error = get_residue_identifiers(structure, model_number=model_number, author_chain=autha['chain'], author_resnum=int(autha['number']), author_resname=autha['name'], author_ins_code=insa)
    if error:
      print(error, file=sys.stderr)
      continue

    resultb, error = get_residue_identifiers(structure, model_number=model_number, author_chain=authb['chain'], author_resnum=int(authb['number']), author_resname=authb['name'], author_ins_code=insb)
    if error:
      print(error, file=sys.stderr)
      continue

    # seq_id can be None for pairs to non-polymer nucleotides
    if not isinstance(resulta['seq_id'], int) or not isinstance(resultb['seq_id'], int):
      print(f"Skipping pair with invalid seq_id: {auth_key(auth1)} {auth_key(auth2)} (seq_id_1={resulta['seq_id']}, seq_id_2={resultb['seq_id']})", file=sys.stderr)
      continue

    # add both directions to processed
    processed.append(key12)
    processed.append(auth_key(auth2) + "_" + auth_key(auth1))

    base_pair_id += 1
    annotation_id += 1

    # populate _ndb_base_pair_list loop (struct_oper_id is always identity, like in original script
    ndb_base_pair_list_loop.add_row([str(base_pair_id), str(model_number), resulta['asym_id'], resulta['entity_id'], str(resulta['seq_id']), resulta['comp_id'], resulta['PDB_ins_code'], resulta['alt_id'], identity_id, resultb['asym_id'], resultb['entity_id'], str(resultb['seq_id']), resultb['comp_id'], resultb['PDB_ins_code'], resultb['alt_id'], identity_id, autha['chain'], str(autha['number']), authb['chain'], str(authb['number'])])

    # populate _ndb_base_pair_annotation loop
    ndb_base_pair_annotation_loop.add_row([str(annotation_id), str(base_pair_id), long_names[fam[0]], long_names[fam[1]], long_names[fam[2]], bp_families[fam], fam, fam + '_' + autha['name'] + '-' + authb['name'], fam + '_' + autha['name'] + '-' + authb['name'] + '_1'])
except Exception as e:
  print(f"error {e=}, {type(e)=}", file=sys.stderr)
  raise

options = gemmi.cif.WriteOptions()
options.align_pairs = 33
options.align_loops = 30

# and save the final mmcif
doc.write_file(out_file, options)

# get also as string just to visualize if they are no errors
cif_in_string = doc.as_string(options)
print(cif_in_string)
