"""Append an annotator to the _software category of a base-pairing mmCIF."""
import gemmi

# columns of a new _software loop (used when the structure has no _software category)
NEW_COLUMNS = ['pdbx_ordinal', 'name', 'version', 'type', 'location', 'description']


def _cell(value):
  return '?' if value is None else gemmi.cif.quote(value)


def add_software_row(block, software, description):
  """Add one row for `software` (dict with name, commit, repository) to the _software category.

  An existing loop is extended: columns without a value get '?' and pdbx_ordinal is the next
  free integer. Without the category a new loop is created. A category written as key-value
  pairs cannot be extended and is left unchanged. Returns True when a row was added.
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
    return True
  if table.loop is None:
    return False
  ordinals = [int(v) for v in block.find_values('_software.pdbx_ordinal') if v.isdigit()]
  values['pdbx_ordinal'] = str(max(ordinals, default=0) + 1)
  table.append_row([_cell(values.get(tag.split('.', 1)[1].lower())) for tag in table.tags])
  return True
