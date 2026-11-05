"""Item names and class labels of the NDB base-pair categories, as defined in docs/ndb-bp-ext.dic.

LIST_ITEMS and ANNOTATION_ITEMS are the items the exporters write, in the order of the dictionary.
pair_class() builds the `class` label of a pair from the parent bases of its two residues: a modified
residue keeps its own identity in _ndb_base_pair_list and is labelled with its parent base in `class`
(2MG paired with C is cWW_G-C). The class values of the dictionary use the bases A, C, G and U only: a pair
with another parent base (T of DNA, N when unknown) has class '?'.
"""
import re

import gemmi

LIST_ITEMS = [
  'id', 'PDB_model_num',
  'label_asym_id_1', 'label_entity_id_1', 'label_seq_id_1', 'label_comp_id_1',
  'PDB_ins_code_1', 'label_alt_id_1', 'struct_oper_id_1',
  'label_asym_id_2', 'label_entity_id_2', 'label_seq_id_2', 'label_comp_id_2',
  'PDB_ins_code_2', 'label_alt_id_2', 'struct_oper_id_2',
  'auth_asym_id_1', 'auth_seq_id_1', 'auth_asym_id_2', 'auth_seq_id_2',
]

ANNOTATION_ITEMS = ['id', 'base_pair_id', 'orientation', 'base_1_edge', 'base_2_edge',
                    'l-w_family_num', 'l-w_family_name', 'class', 'sub_class']

# names used before docs/ndb-bp-ext.dic: still read by parse_basepairs, never written
LEGACY_LIST_ITEMS = {'id': 'base_pair_id', 'PDB_model_num': 'PDB_model_number'}
for _side in ('1', '2'):
  for _new, _old in (('label_asym_id', 'asym_id'), ('label_entity_id', 'entity_id'), ('label_seq_id', 'seq_id'),
                     ('label_comp_id', 'comp_id'), ('label_alt_id', 'alt_id')):
    LEGACY_LIST_ITEMS[f'{_new}_{_side}'] = f'{_old}_{_side}'
LEGACY_ANNOTATION_ITEMS = {'l-w_family_name': 'l-w_family', 'sub_class': 'subclass'}

STANDARD = {'A': 'A', 'C': 'C', 'G': 'G', 'U': 'U', 'T': 'T', 'DA': 'A', 'DC': 'C', 'DG': 'G', 'DT': 'T', 'DU': 'U'}
BASE_ORDER = 'AGCU'                  # the bases of the class values, in the order of the dictionary (A-G, G-C, C-U, ...)
NUCLEOTIDE_PHOSPHATE = re.compile(r'^([ACGUT])[MDT]P$')    # AMP, GDP, GTP, ...: a 5'-terminal residue of a transcript
SYMMETRIC_FAMILIES = {'cWW', 'tWW', 'cHH', 'tHH', 'cSS', 'tSS'}      # both edges alike: the class lists the bases in BASE_ORDER


def _category(block, name):
  return {tag.lower(): column for tag, column in block.get_mmcif_category(name).items()}


def parent_base_lookup(block):
  """Function comp_id -> one-letter parent base (A, C, G, U, T; N when unknown).

  Order of the sources: the standard residue names; _pdbx_struct_mod_residue.parent_comp_id and
  _chem_comp.mon_nstd_parent_comp_id of the structure; the nucleotide table of gemmi (the letter of the
  parent base of a modified nucleotide, e.g. PSU -> U); the nucleoside phosphates (GTP -> G).
  """
  parents = {}
  for category, comp_item, parent_item in (('_pdbx_struct_mod_residue.', 'label_comp_id', 'parent_comp_id'),
                                           ('_chem_comp.', 'id', 'mon_nstd_parent_comp_id')):
    table = _category(block, category)
    for comp, parent in zip(table.get(comp_item, []), table.get(parent_item, [])):
      if isinstance(parent, str) and parent:
        parents.setdefault(comp, parent)

  def letter(comp_id):
    if comp_id in STANDARD:
      return STANDARD[comp_id]
    if parents.get(comp_id) in STANDARD:
      return STANDARD[parents[comp_id]]
    info = gemmi.find_tabulated_residue(comp_id)
    if info is not None and info.is_nucleic_acid():
      return info.one_letter_code.strip().upper() or 'N'
    phosphate = NUCLEOTIDE_PHOSPHATE.match(comp_id)
    return phosphate.group(1) if phosphate else 'N'

  return letter


def pair_class(family, comp_id_1, comp_id_2, letter):
  """class label of a pair, e.g. cWW_G-C, or '?' when a parent base is not in BASE_ORDER.

  `letter` is the function from parent_base_lookup.
  """
  first, second = letter(comp_id_1), letter(comp_id_2)
  if first not in BASE_ORDER or second not in BASE_ORDER:
    return '?'
  if family in SYMMETRIC_FAMILIES and _rank(first) > _rank(second):
    first, second = second, first
  return f'{family}_{first}-{second}'


def _rank(base):
  return BASE_ORDER.index(base)
