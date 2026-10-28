# FR3D and RNApolis on the six stored examples

Pair-level comparison of the two annotators on `data/tests/Inputs/{8bwt,8vjt,8xzn,9e5i,9hrf,9sfq}.cif`,
produced by the pipeline (T1 MAXIT, T2 FR3D and RNApolis, T3 validation).

## Reproduce

```bash
# from workstreams/ws1-annotation-validation, with the submodules initialised
for id in 8bwt 8vjt 8xzn 9e5i 9hrf 9sfq; do
  nextflow run main.nf --input ../../data/tests/Inputs/$id.cif --outdir out/$id --visualize
done
python bin/comparison_table.py out/8bwt out/8vjt out/8xzn out/9e5i out/9hrf out/9sfq
```

## Result

Pairs are matched on model, both residues (label chain, label sequence number, insertion code) and
symmetry operator; the model compared is the lowest one present in both tools (model 1). Column
definitions are in the "Validation report" section of the README.

| case | model | fr3d pairs | rnapolis pairs | fr3d skipped | rnapolis skipped | distinct | shared | agreed | discordant | fr3d only | rnapolis only | name mismatch |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 8bwt | 1 | 12 | 12 | 0 | 0 | 12 | 12 | 12 | 0 | 0 | 0 | 0 |
| 8vjt | 1 | 13 | 1 | 36 | 0 | 14 | 0 | 0 | 0 | 13 | 1 | 0 |
| 8xzn | 1 | 23 | 25 | 0 | 0 | 27 | 21 | 21 | 0 | 2 | 4 | 0 |
| 9e5i | 1 | 54 | 70 | 0 | 0 | 72 | 52 | 48 | 4 | 2 | 18 | 0 |
| 9hrf | 1 | 32 | 33 | 0 | 2 | 33 | 32 | 32 | 0 | 0 | 1 | 0 |
| 9sfq | 1 | 35 | 40 | 0 | 0 | 40 | 35 | 33 | 2 | 0 | 5 | 0 |

"agreed" counts pairs with the same Leontis-Westhof family in both tools. The comparison measures
agreement between the two annotators, not accuracy against a reference.

## Notes per case

- **8bwt:** a structure with 10 models. FR3D annotates all of them (120 pairs), RNApolis the first
  model only (12 pairs). Only model 1 is compared; the 108 FR3D rows of the other models are reported as
  `rows_outside_compared_model` in `validation_report.json`.
- **8vjt:** all 17 FR3D pairs involve a symmetry operator (contacts between symmetry mates), while
  RNApolis annotates the asymmetric unit only and reports a single pair, so the two sets do not
  overlap. 36 of the 70 FR3D input rows are skipped because the two residues have different
  alternate positions (`alt_position_mismatch`). The 17 written FR3D rows include four contacts
  listed once per conformer (A/A and B/B); the alternate position is not part of the key, so each
  counts once (`duplicate_keys`) and 13 distinct pairs remain.
- **9hrf:** two RNApolis pairs are skipped because a residue has no `label_seq_id`
  (`no_label_seq_id`).
- **9e5i, 9sfq:** the discordant pairs have the same residues and different families.

## Reproduction of the stored files

The full pipeline was run on all six inputs. For every case the regenerated FR3D table, the
regenerated RNApolis JSON (pairs and families) and the `_ndb_base_pair_list` rows of both
base-pair mmCIFs are identical to the files in `data/tests/fr3d_outputs`,
`data/tests/rnapolis_outputs` and `data/tests/conversion_outputs`.

## Versions

| component | version |
|---|---|
| FR3D (`fr3d-python`) | commit `995de6cf5fb67e689ab8b6626988799ec4362f1e`, `pyproject` 0.0.2 |
| RNApolis (`rnapolis-py`) | commit `5ba519b5af6de0086d02e61640bdff56e2ef6bc5`, `pyproject` 0.16.7 |
| MAXIT | 11.400 (bioconda) |
| python | 3.12.16 |
| gemmi | 0.7.5 |
| numpy, scipy, pandas, pulp | 2.5.4, 1.18.1, 3.1.0, 3.3.2 |
| forgi | 2.2.3 |
| Nextflow | 26.04.2 |

`run_manifest.json` of each run records the resolved version of every package.
