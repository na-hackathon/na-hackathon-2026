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
# 1. dev toolchain (Nextflow + Java + gemmi + nf-test + pytest)
mamba env create -f environment.yml      # or: conda env create -f environment.yml
conda activate ws1-dev

# 2. fetch the pinned annotator sources (see "Annotator submodules" below)
git submodule update --init \
    annotations_tools/fr3d-python \
    annotations_tools/rnapolis-py

# 3. smoke test — runs end-to-end on the bundled 9CFN data
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
  validation/validation_report.{json,tsv}          # T3 cross-tool comparison
  validation/<tool>.{csv,tsv,bg,jpg}               # T3 optional, with --visualize
```

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

# Optional T3 export: per-tool CSV/TSV pair tables + forgi BulgeGraph (.bg) + matplotlib .jpg
nextflow run main.nf --input <file> --visualize

# Custom output dir
nextflow run main.nf --input <file> --outdir out/8xzn
```

### Parameters

| Param | Default | Meaning |
|---|---|---|
| `--input`      | _(required)_      | Structure file, `.pdb` or `.cif` |
| `--converter`  | `maxit`           | T1 converter: `maxit` (bioconda) \| `maxit-docker` (tzok/maxit) \| `gemmi` (deprecated) |
| `--annotators` | `fr3d,rnapolis`   | Comma-separated annotator names (one `envs/<tool>.yml` each) |
| `--visualize`  | `false`           | Also export per-tool CSV/TSV/forgi BulgeGraph (`.bg`) + render (`.jpg`) |
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
  nextflow.config          base nf-test config (conda + docker off)
  docker.config            per-test overlay that turns docker.enabled on
  parse.nf.test            PARSE process, 6 PDBs × 2 tools
  annotate.nf.test         ANNOTATE process, 6 PDBs × 2 tools + unsupported-tool guard
  convert.nf.test          CONVERT process, MaxIT on each sample structure
  pipeline.nf.test         full main.nf end-to-end, per input × per converter mode
  test_parse_conversion.py pytest unit test for the bin/ parser scripts
```

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
| `tests/parse.nf.test`         | `gemmi` on PATH (host or `ws1-dev` env) |
| `tests/annotate.nf.test`      | `rnapolis` + `fr3d` deps (`mmcif-pdbx`, etc.); submodules initialized |
| `tests/convert.nf.test`       | `maxit` (bioconda) or the docker daemon for the `tzok/maxit` image |
| `tests/pipeline.nf.test`      | Same as the per-stage suites combined; pulls `tzok/maxit` for docker tests |
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

A stage is wired in `modules/*.nf` and calls a CLI on `PATH` (Nextflow puts `bin/` there).
Teams replace the CLI body; **keep the signature stable**.

| Stage | CLI | Signature |
|---|---|---|
| T1 Convert  | `bin/ws1-convert`  | `--input <pdb\|cif> --out <cif> [--converter maxit\|maxit-docker\|gemmi]` |
| T2 Annotate | `bin/ws1-annotate` | `--tool <name> --structure <cif> --out <prefix> --src <annotations_tools dir>` |
| T2 Parse    | `bin/ws1-parse`    | `--tool <name> --raw <file> --structure <cif> --out <cif>` |
| T3 Validate | `bin/ws1-validate` | `--annotations <cif…> --out <prefix>` |
| T3 Visualize | `bin/export_basepairs.py` | `<basepairs.cif> <out_prefix>` — emits `.csv`, `.tsv`, `.bg`, `.jpg` |

## Conversion (T1)

`modules/convert.nf` dispatches via `--converter`:

- **`maxit`** (default) — MAXIT from `bioconda::maxit`, installed into `envs/convert.yml`.
  PDB inputs run as `maxit -o 1`, mmCIF inputs as `maxit -o 8` (refresh).
- **`maxit-docker`** — same MAXIT, but from the prebuilt
  [`tzok/maxit`](https://github.com/tzok/maxit-docker) image. Useful as a fallback when
  bioconda::maxit can't be installed (e.g. ARM Macs).
- **`gemmi`** (DEPRECATED) — lightweight conda fallback for development only.

## Adding an annotator

1. Add the tool's source to `annotations_tools/` (preferably as a submodule pin).
2. Add `envs/<tool>.yml` with its dependencies.
3. Add a dispatch branch for it in `bin/ws1-annotate` (FR3D and RNApolis are the templates).
4. Teach `bin/ws1-parse` to read its native output (drop a small `<tool>_to_*` script in `bin/`).
5. Include the name in `--annotators` — the `ANNOTATE` process is generic over the tool.

## Layout

```
main.nf            workflow wiring (Input → T1 → T2 → T3)
nextflow.config    params + conda/docker enabled by default
modules/           one .nf process per stage
bin/               stage CLIs (ws1-*) + the python parser scripts they call
envs/              per-process conda envs (convert/fr3d/rnapolis/parse/visualize)
annotations_tools/ pinned annotator submodules (fr3d-python, rnapolis-py)
tests/             nf-test + pytest suites (see Tests above)
environment.yml    ws1-dev developer toolchain (NOT a pipeline env)
```

## Status

End-to-end runnable: every stage has a real implementation (no stubs left), all four
nf-test suites + the pytest suite are green on the bundled test inputs, and both
`--converter maxit` (conda) and `--converter maxit-docker` paths are exercised by
`tests/pipeline.nf.test`.
