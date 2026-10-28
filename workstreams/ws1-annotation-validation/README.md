# Workstream 1 — Annotation & validation workflow

A reproducible **Nextflow** pipeline that takes a nucleic-acid structure, runs base-pair
annotators, standardizes their output to base-pairing mmCIF, and validates / compares the
results across tools.

```
Input ──▶ T1 Conversion ──▶ T2 Core (annotate ▸ parse) ──▶ T3 Exploration (validate)
          (any → mmCIF)      FR3D · RNApolis · …            compare across tools
```

## Quickstart

```bash
# 1. toolchain (Nextflow, Java, gemmi, nf-test, pytest)
mamba env create -f environment.yml      # or: conda env create -f environment.yml
conda activate ws1-dev

# 2. fetch the pinned annotator sources (see "Annotator submodules" below)
git submodule update --init \
    annotations_tools/fr3d-python \
    annotations_tools/rnapolis-py

# 3. smoke test — runs end-to-end on the bundled 9CFN data
#    (the first run builds the pipeline env envs/ws1.yml, about 2.7 GB; later runs reuse it)
nextflow run main.nf -profile test

# 4. a real run
nextflow run main.nf --input ../../data/tests/Inputs/8xzn.cif
```

### Annotator submodules

FR3D and RNApolis are pinned as git submodules under `annotations_tools/`. `ws1-annotate`
runs them straight from those checkouts (no `pip install`) via `PYTHONPATH`, so the exact
commit recorded in the parent repo is what executes.

```bash
# from the repo root (na-hackathon-2026/), once after cloning
git submodule update --init \
    workstreams/ws1-annotation-validation/annotations_tools/fr3d-python \
    workstreams/ws1-annotation-validation/annotations_tools/rnapolis-py

# or, equivalently, from this workstream's directory
git submodule update --init \
    annotations_tools/fr3d-python \
    annotations_tools/rnapolis-py

# alternative: clone the parent repo with all submodules in one shot
git clone --recurse-submodules <repo-url>

# in an existing clone that didn't recurse:
git submodule update --init --recursive

# later, pull the pinned commits if upstream bumps them
git submodule update --remote \
    annotations_tools/fr3d-python \
    annotations_tools/rnapolis-py
```

A correctly initialised tree should have these files present:

```
annotations_tools/fr3d-python/fr3d/classifiers/NA_pairwise_interactions.py
annotations_tools/rnapolis-py/src/rnapolis/annotator.py
```

If either is missing, ANNOTATE will fail with `python: can't open file ...` — run the
`git submodule update --init` command above.

Results are written under `--outdir` (default `results/`):

```
results/
  mmcif/            <name>.std.cif                 # T1 standardized structure
  annotations/<tool>/<tool>_raw.*                  # T2 raw tool output
  basepairs/<tool>_basepairs.cif                   # T2 standardized base pairs
  basepairs/<tool>_basepairs.skipped.tsv           # T2 input pairs not written, with the reason
  validation/validation_report.{json,tsv}          # T3 cross-tool comparison
  validation/<tool>.{csv,tsv,bg,jpg}               # T3 optional, with --visualize
```

### Validation report

`validation/validation_report.{json,tsv}` compares the tools pair by pair when more than one
annotator ran. Two tools report the same pair when the model, both residues (label chain, label
sequence number, insertion code) and the symmetry operator match; residue names are compared
separately.

| Field | Meaning |
|---|---|
| `rows` | pairs in each tool's `<tool>_basepairs.cif` |
| `pair_counts` | distinct pairs per tool, after the model selection |
| `total_distinct_pairs` | pairs found by at least one tool |
| `shared_by_all_tools` | pairs found by every tool |
| `agreed_by_all_tools` | shared pairs with the same Leontis–Westhof family in every tool |
| `discordant_family` | shared pairs whose family differs between tools |
| `only_in` | pairs found by a single tool, per tool |
| `comp_id_mismatch` | shared pairs whose residue names differ between tools |
| `duplicate_keys` | rows of a tool that fall on an already counted pair (e.g. alternate positions) |
| `model`, `rows_outside_compared_model` | model compared, and rows of other models left out |
| `rows_with_symmetry_operator` | rows involving a symmetry-related residue |

By default the lowest model present in every tool is compared (RNApolis annotates the first
model only). Other matching rules are available by running the script directly:

```bash
python bin/compare_basepairs.py --compare --out report \
    --annotations results/basepairs/fr3d_basepairs.cif results/basepairs/rnapolis_basepairs.cif \
    [--model all|N] [--asu-only] [--key auth]
```

`--model all|N` selects the models, `--asu-only` drops pairs that involve a symmetry mate, and
`--key auth` matches on auth chain, auth number and residue name only (no model, insertion code
or symmetry operator).

### Skipped pairs

`basepairs/<tool>_basepairs.skipped.tsv` lists the annotator rows that were not written to the
base-pair mmCIF, with columns `reason`, the two residues, the family and `detail`. Reasons:
`unsupported_family` (for example FR3D "near" pairs), `alt_position_mismatch`,
`symmetry_code_not_found`, `residue_not_found` and `no_label_seq_id`. A pair listed in both
directions, or already written, is counted but not listed. Each exporter prints a summary line:
rows read = pairs written + skipped + duplicates.

## Pipeline commands

The pipeline accepts an input structure positionally or via `--input`, lets you pick which
annotators to run, and which MAXIT path to use for T1 conversion.

```bash
# Minimal (default: maxit via bioconda, both annotators, no extra exports)
nextflow run main.nf --input <structure.pdb|cif>

# Single annotator
nextflow run main.nf --input <file> --annotators rnapolis
nextflow run main.nf --input <file> --annotators fr3d

# Both annotators -> VALIDATE step also runs
nextflow run main.nf --input <file> --annotators fr3d,rnapolis

# Pick the converter path
nextflow run main.nf --input <file> --converter maxit          # default, bioconda maxit
nextflow run main.nf --input <file> --converter maxit-docker   # tzok/maxit container
nextflow run main.nf --input <file> --converter gemmi          # DEPRECATED fallback

# Optional T3 export: per-tool CSV/TSV pair tables, and per nucleic-acid chain a forgi
# BulgeGraph (.bg) and a matplotlib render (.jpg)
nextflow run main.nf --input <file> --visualize

# Custom output dir
nextflow run main.nf --input <file> --outdir out/8xzn
```

### Parameters

| Param | Default | Meaning |
|---|---|---|
| `--input`      | _(required)_      | Structure file, `.pdb` or `.cif` |
| `--converter`  | `maxit`           | T1 converter: `maxit` (bioconda) \| `maxit-docker` (tzok/maxit) \| `gemmi` (deprecated) |
| `--annotators` | `fr3d,rnapolis`   | Comma-separated annotator names (those wired in `bin/ws1-annotate`) |
| `--visualize`  | `false`           | Also export per-tool CSV/TSV tables and, per chain, a forgi BulgeGraph (`.bg`) and render (`.jpg`) |
| `--outdir`     | `results`         | Output directory |

`-profile test` sets `--input` to the bundled `data/9CFN/9CFN_clean.cif`. Conda and Docker
are both turned on globally in `nextflow.config`; pick the converter, the rest follows.

## Tests

Two test stacks live under `tests/`:

- **nf-test** — drives the Nextflow modules / pipeline (one `.nf.test` per stage + one for
  the full pipeline).
- **pytest** — drives the bin/ Python parsers (`fr3d_basepair_to_json.py`,
  `json_fr3d_mmcif.py`, `json_rnapolis_mmcif.py`) against reference outputs.

```
tests/
  nextflow.config          base nf-test config (conda on, env cache shared across tests; docker off)
  docker.config            per-test overlay that turns docker.enabled on
  parse.nf.test            PARSE process, 6 PDBs × 2 tools
  annotate.nf.test         ANNOTATE process, 6 PDBs × 2 tools + unsupported-tool guard
  convert.nf.test          CONVERT process, MaxIT on each sample structure
  pipeline.nf.test         full main.nf end-to-end, per input × per converter mode
  test_parse_conversion.py pytest unit test for the bin/ parser scripts
```

Run the tests from an activated `ws1-dev` (`conda activate ws1-dev`). nf-test builds the
pipeline env `envs/ws1.yml` into `.nf-test/conda` during the first test and reuses it afterwards.

If the annotator tests fail with `No module named 'numpy'`, `ws1-dev/bin` precedes the process
env on `PATH`: activate `ws1-dev` instead of prepending its `bin/` to `PATH`.

### Run individual nf-test suites

```bash
nf-test test tests/parse.nf.test          # T2 parsers (12 tests)
nf-test test tests/annotate.nf.test       # T2 annotators (13 tests)
nf-test test tests/convert.nf.test        # T1 conversion (6 tests)
nf-test test tests/pipeline.nf.test       # end-to-end main.nf (15+ tests)
```

### Run everything

```bash
nf-test test tests/         # all nf-tests
pytest tests/               # all pytests
```

### What each suite needs

| Suite | Needs |
|---|---|
| `tests/parse.nf.test`         | the `ws1-dev` toolchain (the process env `envs/ws1.yml` is built on first use) |
| `tests/annotate.nf.test`      | same; annotator submodules initialized |
| `tests/convert.nf.test`       | same (maxit comes from `envs/ws1.yml`) |
| `tests/pipeline.nf.test`      | same; the docker tests also need the docker daemon and pull `tzok/maxit` |
| `tests/test_parse_conversion.py` | `gemmi` + `pytest` |

The pipeline tests are parameterized over every structure under `data/tests/Inputs/` (six
mmCIF inputs) and `data/tests/Inputs/pdb/` (the matching PDB inputs), each run twice — once
with `--converter maxit` (conda) and once with `--converter maxit-docker` (tzok/maxit image).

```bash
# Only the conda-side pipeline tests (no docker daemon required)
nf-test test tests/pipeline.nf.test --tag '*conda*'

# Only the docker-side pipeline tests (needs docker access)
nf-test test tests/pipeline.nf.test --tag '*docker*'
```

## The module contract

Each stage is a process in `modules/*.nf` that calls a CLI on `PATH` (Nextflow adds `bin/`).
The CLI signatures below are the interface between the processes and the implementations;
keep them unchanged when modifying a stage.

| Stage | CLI | Signature |
|---|---|---|
| T1 Convert  | `bin/ws1-convert`  | `--input <pdb\|cif> --out <cif> [--converter maxit\|maxit-docker\|gemmi]` |
| T2 Annotate | `bin/ws1-annotate` | `--tool <name> --structure <cif> --out <prefix> --src <annotations_tools dir>` |
| T2 Parse    | `bin/ws1-parse`    | `--tool <name> --raw <file> --structure <cif> --out <cif>` |
| T3 Validate | `bin/ws1-validate` | `--annotations <cif…> --out <prefix>` |
| T3 Visualize | `bin/export_basepairs.py` | `<basepairs.cif> <out_prefix>` — emits `.csv`, `.tsv`, `.bg`, `.jpg` |

## Conversion (T1)

`modules/convert.nf` dispatches via `--converter`:

- **`maxit`** (default) — MAXIT from `bioconda::maxit`, installed into `envs/ws1.yml`.
  PDB inputs run as `maxit -o 1`, mmCIF inputs as `maxit -o 8` (refresh).
- **`maxit-docker`** — same MAXIT, but from the prebuilt
  [`tzok/maxit`](https://github.com/tzok/maxit-docker) image. Useful as a fallback when
  bioconda::maxit can't be installed (e.g. ARM Macs).
- **`gemmi`** (DEPRECATED) — lightweight conda fallback for development only.

## Adding an annotator

1. Add the tool's source to `annotations_tools/` (preferably as a submodule pin).
2. Add its dependencies to `envs/ws1.yml` (the single env shared by every process).
3. Add a dispatch branch for it in `bin/ws1-annotate` (FR3D and RNApolis are the templates).
4. Teach `bin/ws1-parse` to read its native output (drop a small `<tool>_to_*` script in `bin/`).
5. Include the name in `--annotators` — the `ANNOTATE` process is generic over the tool.

## Layout

```
main.nf            workflow wiring (Input → T1 → T2 → T3)
nextflow.config    params + conda/docker enabled by default
modules/           one .nf process per stage
bin/               stage CLIs (ws1-*) + the python parser scripts they call
envs/ws1.yml       the single conda env shared by every process (pinned python, pulp)
annotations_tools/ pinned annotator submodules (fr3d-python, rnapolis-py)
tests/             nf-test + pytest suites (see Tests above)
environment.yml    ws1-dev toolchain: Nextflow, nf-test, pytest, gemmi (not used by pipeline processes)
```

## Status

End-to-end runnable: every stage has a real implementation (no stubs left), all four
nf-test suites + the pytest suite are green on the bundled test inputs, and both
`--converter maxit` (conda) and `--converter maxit-docker` paths are exercised by
`tests/pipeline.nf.test`.
