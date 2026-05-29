#!/usr/bin/python3

import csv
import json
import sys

basepairs = {}
basepairs["summary"] = {}
basepairs["details"] = []

try:
  in_file = sys.argv[1]
  out_file = sys.argv[2] if len(sys.argv) > 2 else 'basepairs.json'
  with open(in_file) as csvfile:
    reader = csv.reader(csvfile, delimiter='\t')
    for row in reader:
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
except:
  print("usage:", sys.argv[0], "fr3d_basepairs.tsv [out.json]", file=sys.stderr)
  raise
