#!/usr/bin/python3

import csv
import json
import sys

basepairs = {}
basepairs["summary"] = {}
basepairs["details"] = []

try:
  pdbid=""
  print(pdbid)
  with open(sys.argv[1]) as csvfile:
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

    with open('basepairs.json', 'w') as f:
      json.dump(basepairs, f)
except:
  print("Error processing basepair file", file=sys.stderr)
  raise
