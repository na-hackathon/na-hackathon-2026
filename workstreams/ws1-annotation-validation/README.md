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
#    (the first run builds the pipeline env envs/ws1.yml, about 2.7 GB, into ~/.cache/ws1-conda;
#     all later runs, nf-test and the API reuse it; NXF_CONDA_CACHEDIR changes the location)
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
  annotations/<tool>/<tool>_provenance.json        # T2 annotator commit, command line, input checksum
  basepairs/<tool>_basepairs.cif                   # T2 standardized base pairs
  basepairs/<tool>_basepairs.skipped.tsv           # T2 input pairs not written, with the reason
  validation/validation_report.{json,tsv}          # T3 cross-tool comparison
  validation/<tool>.{csv,tsv,bg,jpg}               # T3 optional, with --visualize
  run_manifest.json                                # pipeline commit, versions, checksums of the run
  pipeline_info/{trace.txt,report.html,timeline.html}   # Nextflow execution reports
```

The end of a run prints a summary of the validation report (tools, pairs per tool, distinct,
agreed, discordant and single-tool pairs, and the path of the full table).

### Provenance

- **`_software` in `basepairs/<tool>_basepairs.cif`:** one row per annotator is appended to the
  structure's `_software` category (name, `version` = git commit of the annotator submodule,
  `type` = `program`, repository URL as `location`, and `pdbx_ordinal` = next free integer). The
  rows of the input structure are kept.
- **`annotations/<tool>/<tool>_provenance.json`:** annotator name, repository, git commit,
  `git describe`, version from its `pyproject.toml`, the command line, name and sha256 of the
  structure it ran on, and the python version. Values that cannot be determined are `null`.
- **`run_manifest.json`:** commit of the pipeline directory and whether it has uncommitted changes,
  sha256 of `envs/ws1.yml`, Nextflow version, command line, run name, profile, start time and
  parameters, name and sha256 of the input, the provenance of every annotator, and the version of
  every conda package and python distribution in the pipeline env.

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

### Visualization (`--visualize`)

For each annotator, `validation/<tool>.csv` and `.tsv` list the pairs; `<tool>.bg` and `<tool>.jpg`
show the cWW secondary structure, one forgi BulgeGraph and one panel per nucleic-acid chain.

- Only the first model and the pairs inside the asymmetric unit are drawn: pairs of other models
  and pairs with a symmetry mate are not.
- Residues are matched on auth number and insertion code. Modified nucleotides are drawn with the
  one-letter code of their parent base (`N` when unknown); ligands are not drawn.
- Crossing cWW pairs (pseudoknots) are written with `[]` and then `<>` brackets. A residue that is in
  several cWW pairs keeps the first one, and pairs beyond three crossing layers are left out; both are
  reported on stderr.
- Pairs between two chains are not drawn: each chain is a separate graph.

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
| `--env_file`   | `envs/ws1.yml`    | Conda env of every process; `envs/ws1.lock.yml` pins the version of every package |

`-profile test` sets `--input` to the bundled `data/9CFN/9CFN_clean.cif`. Every process runs in the
conda env. Docker is optional: only `--converter maxit-docker` (and the nf-tests tagged `docker`) run
a container; `docker.enabled` is set in `nextflow.config` for that, and no process uses a container
otherwise, so the default pipeline runs on a machine without Docker.

### Environment versions

`envs/ws1.yml` pins the direct dependencies to minor versions. `envs/ws1.lock.yml` lists the exact
version of every package, transitive dependencies included (linux-64, conda and pip). To run with
it: `nextflow run main.nf --env_file envs/ws1.lock.yml`. To regenerate it after changing
`envs/ws1.yml`, create the env from `envs/ws1.yml` (`conda env create -p <prefix> -f envs/ws1.yml`)
and run `bin/ws1-lock <prefix> > envs/ws1.lock.yml`. `run_manifest.json` of every run records the
versions installed.

The conda env is built once into `~/.cache/ws1-conda` (or `$NXF_CONDA_CACHEDIR`) and shared by
runs started from different directories, by nf-test and by the API.

### Concurrent runs

- **Build the env first.** `bin/ws1-env` builds the env into the cache (name `ws1-<hash of the env
  file>`) and prints its prefix. Builds are serialised with `flock`, so concurrent calls build it once,
  and an interrupted build is detected and redone. Nextflow uses an existing env directory as it is:
  `nextflow run main.nf --env_file "$(bin/ws1-env)"`. The API calls `bin/ws1-env` before every run.
- **One output directory per run.** Runs that share `--outdir` overwrite each other's results. Runs
  started in the same project directory share `work/` and `.nextflow/`; give each its own `-w <dir>`
  when they overlap, as the API does.
- Without `bin/ws1-env`, Nextflow builds the env itself on the first run and takes a lock; if that
  build is killed, delete the `env-*` directory it was creating in the cache.

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
pipeline env `envs/ws1.yml` into the shared conda cache during the first test and reuses it afterwards.
`bin/ws1-nftest [nf-test test arguments]` runs `nf-test test` with its log in `.nf-test/nf-test.log`;
plain `nf-test test` writes `.nf-test.log` to the project root, as nf-test has no setting for it.

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
bin/ws1-nftest tests/pipeline.nf.test --tag conda

# Only the docker-side pipeline tests (needs docker access)
bin/ws1-nftest tests/pipeline.nf.test --tag docker
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
  bioconda::maxit can't be installed (e.g. ARM Macs). This is the only converter that needs Docker;
  the annotators, parser and validation still run in the conda env.
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
envs/ws1.yml       the single conda env shared by every process (exact versions)
docs/              six_example_comparison.md: FR3D vs RNApolis on the six stored structures
annotations_tools/ pinned annotator submodules (fr3d-python, rnapolis-py)
tests/             nf-test + pytest suites (see Tests above)
environment.yml    ws1-dev toolchain: Nextflow, nf-test, pytest, gemmi (not used by pipeline processes)
```

## Status

Every stage runs end to end. `tests/` has nf-test suites for each stage and for the full
pipeline (`--converter maxit` and `--converter maxit-docker`) and pytest suites for the Python
code. `docs/six_example_comparison.md` documents the pipeline run on the six stored structures.
