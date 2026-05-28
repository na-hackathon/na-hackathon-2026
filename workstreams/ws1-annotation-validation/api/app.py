"""WS1 control-plane API.

A thin FastAPI layer the React UI drives: advertise options, launch/monitor/cancel
Nextflow runs, stream live progress (SSE), and serve published results.

Run it (from the workstream dir, with the ws1-dev env active):
    uvicorn api.app:app --reload --port 8000

The weblog callback URL must match where this server is reachable; override with
the WS1_API_BASE_URL env var if not http://127.0.0.1:8000.
"""
from __future__ import annotations

import asyncio
import json
import os
import re
from pathlib import Path

from fastapi import FastAPI, File, HTTPException, Request, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse

from .models import RunInfo, RunRequest
from .runner import RunManager, UPLOAD_DIR

BASE_URL = os.environ.get("WS1_API_BASE_URL", "http://127.0.0.1:8000")

# Upload hardening: strict filename charset, structure-only extensions, size cap.
SAFE_FILENAME = re.compile(r"\A[A-Za-z0-9._-]{1,128}\Z")
ALLOWED_SUFFIXES = {".pdb", ".cif", ".mmcif"}
MAX_UPLOAD_BYTES = 50 * 1024 * 1024

app = FastAPI(title="WS1 control plane", version="0.1.0")
# Permissive CORS for dev so the React app (different origin) can call the API.
# Tighten allow_origins before any non-local deployment.
app.add_middleware(
    CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"],
)
mgr = RunManager(BASE_URL)


@app.get("/tools")
def tools() -> dict:
    return {"annotators": mgr.list_annotators()}


@app.get("/profiles")
def profiles() -> dict:
    return {"profiles": mgr.list_profiles()}


@app.post("/uploads")
async def upload(file: UploadFile = File(...)) -> dict:
    name = os.path.basename(file.filename or "")
    if (name in {"", ".", ".."} or not SAFE_FILENAME.match(name)
            or Path(name).suffix.lower() not in ALLOWED_SUFFIXES):
        raise HTTPException(400, "filename must be [A-Za-z0-9._-] ending in .pdb/.cif/.mmcif")
    data = await file.read()
    if len(data) > MAX_UPLOAD_BYTES:
        raise HTTPException(413, "uploaded file too large")
    dest = (UPLOAD_DIR / name).resolve()
    if dest.parent != UPLOAD_DIR.resolve():        # must land directly in uploads
        raise HTTPException(400, "invalid filename")
    dest.write_bytes(data)
    return {"input": str(dest)}


@app.post("/runs", response_model=RunInfo)
async def create_run(req: RunRequest) -> RunInfo:
    try:
        return await mgr.launch(req)
    except ValueError as exc:
        raise HTTPException(400, str(exc))


@app.get("/runs", response_model=list[RunInfo])
def list_runs() -> list[RunInfo]:
    return mgr.list()


@app.get("/runs/{run_id}", response_model=RunInfo)
def get_run(run_id: str) -> RunInfo:
    info = mgr.get(run_id)
    if info is None:
        raise HTTPException(404, "run not found")
    return info


@app.delete("/runs/{run_id}", response_model=RunInfo)
async def cancel_run(run_id: str) -> RunInfo:
    info = await mgr.cancel(run_id)
    if info is None:
        raise HTTPException(404, "run not found")
    return info


@app.post("/runs/{run_id}/_weblog")
async def weblog(run_id: str, request: Request) -> dict:
    """Internal: Nextflow `-with-weblog` posts run/task events here."""
    await mgr.ingest_weblog(run_id, await request.json())
    return {"ok": True}


@app.get("/runs/{run_id}/events")
async def events(run_id: str) -> StreamingResponse:
    """Server-Sent Events: streams live run/task progress until the run ends."""
    run = mgr.runs.get(run_id)
    if run is None:
        raise HTTPException(404, "run not found")

    async def stream():
        queue: asyncio.Queue[dict] = asyncio.Queue()
        run.subscribers.add(queue)            # register before checking, to avoid a race
        try:
            if run.info.finished_at is not None:
                terminal = {"event": "terminated", "status": run.info.status.value,
                            "exit_code": run.info.exit_code}
                yield f"data: {json.dumps(terminal)}\n\n"
                return
            while True:
                event = await queue.get()
                yield f"data: {json.dumps(event)}\n\n"
                if event.get("event") == "terminated":
                    break
        finally:
            run.subscribers.discard(queue)

    return StreamingResponse(stream(), media_type="text/event-stream")


@app.get("/runs/{run_id}/results")
def results(run_id: str) -> dict:
    if mgr.get(run_id) is None:
        raise HTTPException(404, "run not found")
    return {"files": mgr.results(run_id)}


@app.get("/runs/{run_id}/results/{rel:path}")
def result_file(run_id: str, rel: str) -> FileResponse:
    if mgr.get(run_id) is None:
        raise HTTPException(404, "run not found")
    try:
        path = mgr.result_path(run_id, rel)
    except ValueError:
        raise HTTPException(400, "invalid path")
    if path is None:
        raise HTTPException(404, "file not found")
    return FileResponse(path)
