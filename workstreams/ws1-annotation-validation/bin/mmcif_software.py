"""Append an annotator to the _software and _ndb_base_pair_provenance categories of a base-pairing mmCIF."""
import gemmi

# columns of a new _software loop (used when the structure has no _software category)
NEW_COLUMNS = ['pdbx_ordinal', 'name', 'version', 'type', 'location', 'description']
# _ndb_base_pair_provenance of docs/ndb-bp-ext.dic; software_id points to _software.pdbx_ordinal
PROVENANCE_COLUMNS = ['entry_id', 'software_id', 'software_name', 'software_version', 'software_url',
                      'bp_params_version', 'details']


def _cell(value):
  return '?' if value is None else gemmi.cif.quote(value)


def add_software_row(block, software, description):
  """Add one row for `software` (dict with name, commit, repository) to the _software category.

  An existing loop is extended: columns without a value get '?' and pdbx_ordinal is the next
  free integer. Without the category a new loop is created. A category written as key-value
  pairs cannot be extended and is left unchanged. Returns the pdbx_ordinal of the new row, or None
  when no row was added.
  """
  values = {
    'name': software['name'],
    'version': software.get('commit'),
    'type': 'program',
    'location': software.get('repository'),
    'description': description,
  }
  table = block.find_mmcif_category('_software.')
  if not table:
    block.init_mmcif_loop('_software.', NEW_COLUMNS).add_row(
      [_cell({**values, 'pdbx_ordinal': '1'}[column]) for column in NEW_COLUMNS])
    return 1
  if table.loop is None:
    return None
  ordinals = [int(v) for v in block.find_values('_software.pdbx_ordinal') if v.isdigit()]
  values['pdbx_ordinal'] = str(max(ordinals, default=0) + 1)
  table.append_row([_cell(values.get(tag.split('.', 1)[1].lower())) for tag in table.tags])
  return int(values['pdbx_ordinal'])


def provenance_details(software):
  """details of the _ndb_base_pair_provenance row: the annotator, and its checkout when it has local changes."""
  text = f"base pairs annotated by {software['name']} in the WS1 pipeline"
  if software.get('state') == 'modified':
    text += "; the annotator checkout has uncommitted changes"
  return text


def add_provenance_row(block, software, software_id, details):
  """Add a row to _ndb_base_pair_provenance for the annotator (software_id = its _software.pdbx_ordinal).

  The entry id is _entry.id of the structure (the data block name without it). An existing loop is
  extended; a category written as key-value pairs is left unchanged. Returns True when a row was added.
  """
  entry = block.find_value('_entry.id')
  values = {
    'entry_id': gemmi.cif.as_string(entry) if entry else block.name,
    'software_id': str(software_id),
    'software_name': software['name'],
    'software_version': software.get('commit'),
    'software_url': software.get('repository'),
    'bp_params_version': None,
    'details': details,
  }
  table = block.find_mmcif_category('_ndb_base_pair_provenance.')
  if not table:
    block.init_mmcif_loop('_ndb_base_pair_provenance.', PROVENANCE_COLUMNS).add_row(
      [_cell(values[column]) for column in PROVENANCE_COLUMNS])
    return True
  if table.loop is None:
    return False
  table.append_row([_cell(values.get(tag.split('.', 1)[1])) for tag in table.tags])
  return True
