# Workstream 1 — Annotation & validation workflow

A reproducible **Nextflow** pipeline that takes a nucleic-acid structure, runs base-pair
annotators, standardizes their output to base-pairing mmCIF, and validates/compares the
results across tools.

This directory is the **skeleton** (the blue box in [`workflow.mmd`](workflow.mmd)):
orchestration + environments + tests. The annotation/conversion/validation *modules* plug
into it through a small CLI contract — each team owns the logic behind one CLI.

```
Input ──▶ T1 Conversion ──▶ T2 Core (annotate ▸ parse) ──▶ T3 Exploration (validate)
          (any → mmCIF)      FR3D · RNApolis · …            compare across tools
```

## Quickstart

```bash
# 1. dev toolchain (Nextflow + Java + gemmi + nf-test + pytest)
mamba env create -f environment.yml      # or: conda env create -f environment.yml
conda activate ws1-dev

# 2. smoke test — runs end-to-end on the bundled 9CFN data, no tool installs needed
nextflow run main.nf -profile test

# 3. a real run
nextflow run main.nf --input ../../data/9CFN_clean.cif --annotators fr3d,rnapolis -profile conda
```

Results are written under `--outdir` (default `results/`):

```
results/
  mmcif/            <name>.std.cif                 # T1 standardized structure
  annotations/<tool>/<tool>_raw.*                  # T2 raw tool output
  basepairs/<tool>_basepairs.cif                   # T2 standardized base pairs
  validation/validation_report.{json,tsv}          # T3 cross-tool comparison
```

## Parameters

| Param | Default | Meaning |
|---|---|---|
| `--input` | _(required)_ | Structure file, `.pdb` or `.cif` |
| `--annotators` | `fr3d,rnapolis` | Comma-separated annotator names (one `envs/<name>.yml` each) |
| `--outdir` | `results` | Output directory |

## Environments (profiles)

Software is attached **per process** and selected with `-profile`. Default is **conda**;
`mamba` is an opt-in speed profile (needs the `mamba` binary); `docker`/`singularity` for
reproducible / HPC runs; `test` runs the stubs with no env management.

| `-profile` | Effect |
|---|---|
| _(none)_ / `conda` | `conda.enabled` — build per-process envs with conda (default) |
| `mamba` | as conda, but solve with `mamba` (faster; requires mamba) |
| `docker` | run each process in its container |
| `singularity` | containers via Singularity/Apptainer (HPC; `autoMounts`) |
| `test` | no env management; bundled 9CFN input — fast skeleton check |

Why per-process envs? FR3D, RNApolis, RNAview etc. have incompatible dependencies, so each
gets its own isolated environment instead of one shared one.

## The module contract

A stage is wired in `modules/*.nf` and calls a CLI on `PATH` (Nextflow puts `bin/` there).
Teams replace the CLI body; **keep the signature stable**. The shipped `bin/*` are runnable
stubs so the skeleton works today.

| Stage | CLI | Signature |
|---|---|---|
| T1 Convert | `bin/ws1-convert` | `--input <pdb\|cif> --out <cif>` |
| T2 Annotate | `bin/ws1-annotate` | `--tool <name> --structure <cif> --out <prefix>` |
| T2 Parse | `bin/ws1-parse` | `--tool <name> --raw <file> --structure <cif> --out <cif>` |
| T3 Validate | `bin/ws1-validate` | `--annotations <cif…> --out <prefix>` |

Reference implementations of the FR3D/RNApolis → base-pairing-mmCIF conversion and the
comparison step already exist in [`read_write_mmcif/`](read_write_mmcif/) (from PR #46) and
can be adapted into `ws1-parse` / `ws1-validate`.

## Adding an annotator

1. Implement `bin/ws1-annotate` behaviour for the tool (or branch on `--tool`).
2. Add `envs/<tool>.yml` with its dependencies.
3. Teach `bin/ws1-parse` to read its native output.
4. Add the name to `--annotators` — the `ANNOTATE` process is generic over the tool.

## Layout

```
main.nf            workflow wiring (Input → T1 → T2 → T3)
nextflow.config    params + profiles (conda/mamba/docker/singularity/test)
modules/           one .nf process per stage
bin/               stage CLIs — the team contract (stubs for now)
envs/              per-process conda envs (one per tool/stage)
api/               control-plane API (FastAPI) the UI drives — see api/README.md
environment.yml    ws1-dev developer toolchain (this is NOT a pipeline env)
read_write_mmcif/  reference parser/converter snippets (PR #46)
```

## Status

The skeleton runs end-to-end with stub stages (produces a valid, empty report). Replacing a
`bin/*` body with real logic is all a module team needs to do — no changes to the wiring.
