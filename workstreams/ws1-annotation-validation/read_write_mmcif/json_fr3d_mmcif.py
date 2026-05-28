#!/usr/bin/python3

import sys, json, re
from purine_pyrimidine import *
import gemmi
from pathlib import Path

struct_opers = []

bp_families = { "cWW": '1', "tWW": '2', "cWH": '3', "tWH": '4', "cWS": '5', "tWS": '6', "cHH": '7', "tHH": '8', "cHS": '9', "tHS": '10', "cSS": '11', "tSS": '12' }
long_names = {"W": "Watson-Crick", "H": "Hoogsteen", "S": "Sugar", "t": "trans", "c": "cis"}

# list of standard residues
#stdres = [ "A", "U", "G", "C", "DA", "DT", "DG", "DC" ]

# definition of the complementary pairs (colored separately if cWW)
complementary = [ ("A", "U"), ("G", "C"), ("DA", "DT"), ("DG", "DC") ]

processed = []
allres = []

# fr3d  unit_id reference
#unit_id_parts = ['pdb', 'model', 'chain', 'component_id', 'component_number', 'atom_name', 'alt_id', 'insertion_code', 'symmetry']
# ['1EHZ|1|A|2MG|10', 'cHS', '1EHZ|1|A|G|45', '0']
# ["138D|1|A|DG|1||||11_555", "cWW", "138D|1|A|DC|10", "0"]

base_pair_id = 0
annotation_id = 0
fr3d_syms = []

pattern = re.compile("[ct][WHS][WHS]")

try:
  with open(sys.argv[1], 'r') as file:
    pairs = json.load(file)

  mmcif_file = sys.argv[2]
  out_file = sys.argv[3] if len(sys.argv) > 3 else Path(mmcif_file).stem + '_basepairs.cif'

  # TODO replace structure with something else so we can drop the double reading
  structure = gemmi.read_structure(mmcif_file)

  # there is some error getting correct document from structure, so we read it twice
  # BUG doc = structure.make_mmcif_document()
  doc = gemmi.cif.read(mmcif_file)
  block = doc.sole_block()

  asym_loop = block.get_mmcif_category('_struct_asym.')
  struct_oper_list_loop = block.get_mmcif_category('_pdbx_struct_oper_list.')
  
  num_sym_oper = 1
  identity_id='1'
  if (struct_oper_list_loop):
    tagpos = {'id':None, 'type':None, 'name':None}
    struct_oper_dict = {}
#    struct_oper_list_loop = block.find_loop_item('_pdbx_struct_oper_list.id').loop
    num_sym_oper = len(struct_oper_list_loop['id'])
    # TODO in the real case, we will compare the trans and rot operations between gemmi and the _pdbx_struct_oper_list
    for idx, sym_id in enumerate(struct_oper_list_loop['id']):
      struct_oper_dict[struct_oper_list_loop['name'][idx]] = sym_id
      # identity is not always first (e.g. 5a39 was creative)
      # as_string removes quotes
      if (gemmi.cif.as_string(struct_oper_list_loop['type'][idx]) == 'identity operation'):
        identity_id = str(sym_id)

  ndb_base_pair_list_loop = block.init_mmcif_loop('_ndb_base_pair_list.', ['base_pair_id', 'PDB_model_number', 'asym_id_1', 'entity_id_1', 'seq_id_1', 'comp_id_1', 'PDB_ins_code_1', 'alt_id_1', 'struct_oper_id_1', 'asym_id_2', 'entity_id_2', 'seq_id_2', 'comp_id_2', 'PDB_ins_code_2', 'alt_id_2', 'struct_oper_id_2', 'auth_asym_id_1', 'auth_seq_id_1', 'auth_asym_id_2', 'auth_seq_id_2'])
  ndb_base_pair_annotation_loop = block.init_mmcif_loop('_ndb_base_pair_annotation.', ['id', 'base_pair_id', 'orientation', 'base_1_edge', 'base_2_edge', 'l-w_family_num', 'l-w_family', 'class', 'subclass'])

except:
  print("usage:",sys.argv[0],"basepairs.json","structure.cif","[out.cif]", file=sys.stderr)
  raise

def get_atom_identifiers(structure, model_number, author_chain, author_resnum, author_resname, author_ins_code, author_alt_id):
  try:
    model = next(m for m in structure if m.num == int(model_number))
  except StopIteration:
    return None, f"Model '{model_number}' not found."

  #print("model", model.name, "found", file=sys.stderr)
  for chain in model:
    if chain.name.strip() == author_chain:
      for residue in chain:
        if (residue.seqid.num == author_resnum and residue.name == author_resname and (residue.seqid.icode if residue.seqid.icode != ' ' else ' ') == author_ins_code):
          for atom in residue:
            if (atom.altloc if atom.altloc != '\x00' else '') == author_alt_id:
              return {
                'model': model.num,
                'asym_id': residue.subchain,
                'entity_id': asym_loop["entity_id"][asym_loop["id"].index(residue.subchain)] ,
                'seq_id': residue.label_seq,
                'comp_id': residue.name,
                'PDB_ins_code': (residue.seqid.icode if residue.seqid.icode != ' ' else '?'),
                'alt_id': (atom.altloc if atom.altloc != '\x00' else '.')
              }, None

          return None, f"Residue found but no atom with alt_id '{author_alt_id}' present."

  return None, f"No match found for model '{model_number}', chain '{author_chain}', residue '{author_resnum}', name '{author_resname}', insertion code '{author_ins_code}', alt_id '{author_alt_id}'."

try:
  for pair in pairs["details"]:
    if pair[0]+"_"+pair[2] not in processed:
      # as long as we are using FR3D, we will skip the "near" pairs or anything unexpected
      if ( not pattern.fullmatch(pair[1]) ):
        print("Skipping unsupported pairing family", pair, file=sys.stderr)
        continue
      # check if the family is in the right order (only the bp_families are allowed)
      if (pair[1] in bp_families):
        resa = pair[0].split("|")
        resb = pair[2].split("|")
        fam = pair[1]
      else:
        # change the order of residues and family name
        resb = pair[0].split("|")
        resa = pair[2].split("|")
        fam = pair[1][0]+pair[1][2]+pair[1][1]

      lena = len(resa)
      alta = resa[6] if lena >= 7 else ''
      insa = resa[7] if (lena >= 8 and resa[7] != '') else ' '
      syma = resa[8] if lena > 8 else identity_id
      if (syma not in fr3d_syms):
        fr3d_syms.append(syma)

      lenb = len(resb)
      altb = resb[6] if lenb >= 7 else ''
      insb = resb[7] if (lenb >= 8 and resb[7] != '') else ' '
      symb = resb[8] if lenb > 8 else identity_id
      if (symb not in fr3d_syms):
        fr3d_syms.append(symb)

      if (alta != '' and altb != '' and alta != altb):
        print("Skipping suspicious alt position combination", alta, altb, "in", pair, file=sys.stderr)
        continue

      # TODO to be corected when gemmi allows access to _pdbx_struct_oper_list.name and/or when our base-pairing code is used instead of fr3d
      if (len(fr3d_syms) > num_sym_oper):
        sys.exit("Number of FR3D symmetry operations is larger than size of _pdbx_struct_oper_list", pair)

      if (syma != identity_id):
        if (syma in struct_oper_dict):
          syma = struct_oper_dict[syma]
        else:
          print("Symmetry code", syma, "not found in _pdbx_struct_oper_list for", pair, file=sys.stderr)
          continue

      if (symb != identity_id):
        if (symb in struct_oper_dict):
          symb = struct_oper_dict[symb]
        else:
          print("Symmetry code", symb, "not found in _pdbx_struct_oper_list for", pair, file=sys.stderr)
          continue

      resulta, error = get_atom_identifiers(structure, model_number=str(resa[1]), author_chain=resa[2], author_resnum=int(resa[4]), author_resname=resa[3], author_ins_code=insa, author_alt_id=alta)
      if error:
        print(error, file=sys.stderr)
        continue

      resultb, error = get_atom_identifiers(structure, model_number=str(resb[1]), author_chain=resb[2], author_resnum=int(resb[4]), author_resname=resb[3], author_ins_code=insb, author_alt_id=altb)
      if error:
        print(error, file=sys.stderr)
        continue

      # Validate that seq_id values are valid integers (it can be None for pairs to non-polymer nucleotides)
      if not isinstance(resulta['seq_id'], int) or not isinstance(resultb['seq_id'], int):
        print(f"Skipping pair with invalid seq_id: {pair} (seq_id_1={resulta['seq_id']}, seq_id_2={resultb['seq_id']})", file=sys.stderr)
        continue

      if resa not in allres:
        allres.append(pair[0])

      if resb not in allres:
        allres.append(pair[2])

      # add both directions to processed
      processed.append(pair[0]+"_"+pair[2])
      processed.append(pair[2]+"_"+pair[0])

      # FR3D sometimes returns the same pairs within ASU and then from different symmetry-generated copy
      # try taking each pair only once
      if (syma != identity_id and symb != identity_id and syma == symb):
        print("Same symmetry pair converted to within ASU", pair, file=sys.stderr)
        syma = identity_id
        symb = identity_id
        # add also ASU variants to the processed list
        del resa[8]
        del resb[8]
        for x in range(7,4,-1):
          if (resa[x] == ''):
            del resa[x]
          else:
            break
        for x in range(7,4,-1):
          if (resb[x] == ''):
            del resb[x]
          else:
            break
        if ('|'.join(resa)+"_"+'|'.join(resb) not in processed):
          processed.append('|'.join(resa)+"_"+'|'.join(resb))
        else:
          print("skipF", pair, file=sys.stderr)
          continue
        if ('|'.join(resb)+"_"+'|'.join(resa) not in processed):
          processed.append('|'.join(resb)+"_"+'|'.join(resa))
        else:
          print("skipR", pair, file=sys.stderr)
          continue

      base_pair_id += 1
      annotation_id += 1

      # populate _ndb_base_pair_list loop
      ndb_base_pair_list_loop.add_row([str(base_pair_id), str(resa[1]), resulta['asym_id'], resulta['entity_id'], str(resulta['seq_id']), resulta['comp_id'], resulta['PDB_ins_code'], resulta['alt_id'], syma, resultb['asym_id'], resultb['entity_id'], str(resultb['seq_id']), resultb['comp_id'], resultb['PDB_ins_code'], resultb['alt_id'], symb, resa[2], resa[4], resb[2], resb[4]])

      # populate _ndb_base_pair_annotation loop; subclass hardcoded to class_1
      ndb_base_pair_annotation_loop.add_row([str(annotation_id), str(base_pair_id), long_names[fam[0]], long_names[fam[1]], long_names[fam[2]], bp_families[fam], fam, fam+'_'+resa[3]+'-'+resb[3], fam+'_'+resa[3]+'-'+resb[3]+'_1'])
    else:
      #print(pair[0]+"_"+pair[2], "already processed", file=sys.stderr)
      pass
except Exception as e:
  print(f"error {e=}, {type(e)=}", file=sys.stderr)
  raise

# set "nice output" writing options for gemmi
options = gemmi.cif.WriteOptions()
options.align_pairs = 33
options.align_loops = 30

# write to mmCIF file
doc.write_file(out_file, options)

# get also as string
cif_in_string = doc.as_string(options)
print(cif_in_string)

### base pairing mmCIF categories used here
# TODO add the provenance

# see ndb-bp-ext.dic
# loop_
# _ndb_base_pair_list.base_pair_id
# _ndb_base_pair_list.PDB_model_number
# _ndb_base_pair_list.asym_id_1
# _ndb_base_pair_list.entity_id_1
# _ndb_base_pair_list.seq_id_1
# _ndb_base_pair_list.comp_id_1
# _ndb_base_pair_list.PDB_ins_code_1
# _ndb_base_pair_list.alt_id_1
# _ndb_base_pair_list.struct_oper_id_1 → mandatory → _pdbx_struct_oper_list.id
# _ndb_base_pair_list.asym_id_2
# _ndb_base_pair_list.entity_id_2
# _ndb_base_pair_list.seq_id_2
# _ndb_base_pair_list.comp_id_2
# _ndb_base_pair_list.PDB_ins_code_2
# _ndb_base_pair_list.alt_id_2
# _ndb_base_pair_list.struct_oper_id_2 → mandatory → _pdbx_struct_oper_list.id
# _ndb_base_pair_list.auth_asym_id_1
# _ndb_base_pair_list.auth_seq_id_1
# _ndb_base_pair_list.auth_asym_id_2
# _ndb_base_pair_list.auth_seq_id_2
# 1 1 1 A 1 C . . 1 1 A 12 G . . 1 P 11 P 22

# loop_
# _ndb_base_pair_annotation.id
# _ndb_base_pair_annotation.base_pair_id
# _ndb_base_pair_annotation.orientation
# _ndb_base_pair_annotation.base_1_edge
# _ndb_base_pair_annotation.base_2_edge
# _ndb_base_pair_annotation.l-w_family_num
# _ndb_base_pair_annotation.l-w_family
# _ndb_base_pair_annotation.class
# _ndb_base_pair_annotation.subclass
# 15 15 trans        Sugar        Sugar 12 tSS tSS_A-G tSS_A-G_1
# 16 16   cis Watson-Crick Watson-Crick  1 cWW cWW_G-C cWW_G-C_1
# 17 17 trans    Hoogsteen        Sugar 10 tHS tHS_A-G tHS_A-G_1
# 18 18   cis Watson-Crick Watson-Crick  1 cWW cWW_G-U cWW_G-U_1
# #
