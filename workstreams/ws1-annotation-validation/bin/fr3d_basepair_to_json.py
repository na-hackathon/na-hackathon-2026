#!/usr/bin/python3
"""Convert an FR3D basepair table (TSV: unit id, family, unit id, ...) to JSON for json_fr3d_mmcif.py.

Usage: fr3d_basepair_to_json.py fr3d_basepairs.tsv [out.json]   (default out: basepairs.json)

The JSON has "summary" (rows per family) and "details" (the rows as read).
"""

import csv
import json
import sys

basepairs = {}
basepairs["summary"] = {}
basepairs["details"] = []

if len(sys.argv) < 2:
  sys.exit(f"usage: {sys.argv[0]} fr3d_basepairs.tsv [out.json]")

in_file = sys.argv[1]
out_file = sys.argv[2] if len(sys.argv) > 2 else 'basepairs.json'
with open(in_file) as csvfile:
  reader = csv.reader(csvfile, delimiter='\t')
  for row in reader:
    if not row:
      continue
    if len(row) < 3:
      sys.exit(f"{in_file}:{reader.line_num}: expected at least 3 tab-separated columns (unit id, family, unit id), got {len(row)}")
    #a=row[0].split('|')
    #b=row[2].split('|')
    LW=row[1]
    if (LW not in basepairs["summary"]):
      basepairs["summary"][LW] = 1
    else:
      basepairs["summary"][LW] += 1
    basepairs["details"].append(row)

  with open(out_file, 'w') as f:
    json.dump(basepairs, f)
