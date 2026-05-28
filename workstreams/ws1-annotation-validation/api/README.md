# WS1 control-plane API

A thin **FastAPI** backend the UI drives to control the Nextflow pipeline: list
options, launch / monitor / cancel runs, stream live progress, and serve results.
It wraps the `nextflow` CLI (the skeleton in the parent directory) — there is no
native Nextflow control API, so this layer provides one.

```
React UI ──HTTP/SSE──▶ this API ──subprocess──▶ nextflow run main.nf …
                          ▲                              │
                          └────── -with-weblog (events) ─┘
```

## Run it

```bash
conda activate ws1-dev
cd workstreams/ws1-annotation-validation
uvicorn api.app:app --reload --port 8000
```

Interactive docs (OpenAPI) at `http://127.0.0.1:8000/docs`. If the server is not
reachable at `http://127.0.0.1:8000`, set `WS1_API_BASE_URL` so Nextflow's weblog
callback points back correctly.

## Endpoints

| Method | Path | Purpose |
|---|---|---|
| GET | `/tools` | annotator names available (from `envs/*.yml`) |
| GET | `/profiles` | Nextflow profiles (`conda`, `mamba`, `docker`, `singularity`, `test`) |
| POST | `/uploads` | upload a `.pdb`/`.cif`/`.mmcif`; returns an opaque `upload_id` |
| POST | `/runs` | launch a run — body `{input?: upload_id, annotators?, profile}` |
| GET | `/runs` | list runs |
| GET | `/runs/{id}` | run status + task progress |
| DELETE | `/runs/{id}` | cancel a run |
| GET | `/runs/{id}/events` | **SSE** live progress until the run ends |
| GET | `/runs/{id}/results` | list published result files |
| GET | `/runs/{id}/results/{path}` | download a result file |

`POST /runs/{id}/_weblog` is internal — Nextflow posts events there.

## Example

```bash
# fast end-to-end run on the bundled 9CFN (no upload, no tool installs)
curl -s -X POST localhost:8000/runs -H 'content-type: application/json' \
     -d '{"profile":"test"}'
# -> {"id":"<run_id>", "status":"running", ...}

# run on your own structure: upload first, then pass the returned id
UID=$(curl -s -X POST localhost:8000/uploads -F 'file=@my.cif' | jq -r .upload_id)
curl -s -X POST localhost:8000/runs -H 'content-type: application/json' \
     -d "{\"profile\":\"conda\",\"input\":\"$UID\"}"

curl -s localhost:8000/runs/<run_id>                       # status
curl -sN localhost:8000/runs/<run_id>/events               # live SSE stream
curl -s localhost:8000/runs/<run_id>/results               # {"files":[...]}
curl -s localhost:8000/runs/<run_id>/results/validation/validation_report.json
```

## Security model

No user-controlled value ever forms a filesystem path:

- **Inputs** are referenced by an opaque `upload_id` (a dict key), and the server
  resolves it to a path it generated; arbitrary server paths can't be requested.
- **Uploads** are stored as `<generated-id><allowlisted-suffix>` — the original
  filename never forms the path.
- **Results** are served only if the requested relative path is in that run's own
  result listing (built from a trusted directory walk).

## Limitations (first prototype)

- Run state + the upload registry are in-memory (lost on restart); no persistence/DB yet.
- CORS is wide open for dev — tighten `allow_origins` before any deployment.
- No auth; intended for local/hackathon use.
- Per-run dirs and uploads live under `.api_runs/` (git-ignored).
