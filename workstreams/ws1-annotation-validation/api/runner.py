"""Run management: launch Nextflow, ingest its weblog events, expose run state.

The pipeline itself is the skeleton in the parent directory (main.nf). Each run
gets an isolated dir under .api_runs/<id>/ (work + results). Live progress comes
from Nextflow's `-with-weblog`, which POSTs JSON events back to this service.
"""
from __future__ import annotations

import asyncio
import time
import uuid
from pathlib import Path

from .models import RunInfo, RunStatus, TaskProgress

PROJECT_DIR = Path(__file__).resolve().parent.parent   # the Nextflow project (main.nf)
RUNS_DIR = PROJECT_DIR / ".api_runs"
UPLOAD_DIR = RUNS_DIR / "_uploads"
RESERVED_ENVS = {"convert", "parse"}                   # stage envs, not annotators
PROFILES = ["conda", "mamba", "docker", "singularity", "test"]


class Run:
    def __init__(self, info: RunInfo):
        self.info = info
        self.proc: asyncio.subprocess.Process | None = None
        self.subscribers: set[asyncio.Queue[dict]] = set()   # one queue per SSE client
        # Trusted run dir, built from the server-generated id (never user input).
        self.dir: Path = RUNS_DIR / info.id

    def publish(self, event: dict) -> None:
        for queue in list(self.subscribers):
            queue.put_nowait(event)


class RunManager:
    def __init__(self, base_url: str):
        self.base_url = base_url.rstrip("/")
        self.runs: dict[str, Run] = {}
        self.uploads: dict[str, Path] = {}      # upload_id -> trusted stored path
        UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

    # ---- discovery (so the UI can build its form) ----
    def list_annotators(self) -> list[str]:
        envs = (PROJECT_DIR / "envs").glob("*.yml")
        return sorted(p.stem for p in envs if p.stem not in RESERVED_ENVS)

    def list_profiles(self) -> list[str]:
        return list(PROFILES)

    def register_upload(self, suffix: str, data: bytes) -> str:
        """Store an uploaded structure under a server-generated id; return the id.
        `suffix` is an allowlisted literal, so the path has no user-controlled part."""
        upload_id = uuid.uuid4().hex[:12]
        dest = UPLOAD_DIR / f"{upload_id}{suffix}"
        dest.write_bytes(data)
        self.uploads[upload_id] = dest
        return upload_id

    # ---- run lifecycle ----
    async def launch(self, req) -> RunInfo:
        # Validate every user-supplied value before it reaches the command line.
        if req.profile not in PROFILES:
            raise ValueError(f"unknown profile {req.profile!r}")
        annotators = req.annotators
        if annotators:
            allowed = set(self.list_annotators())
            unknown = [a for a in annotators if a not in allowed]
            if unknown:
                raise ValueError(f"unknown annotator(s): {unknown}")
        input_path = None
        if req.input:                           # req.input is an upload id (used only as a key)
            stored = self.uploads.get(req.input)
            if stored is None:
                raise ValueError("unknown upload id")
            input_path = str(stored)            # trusted path from the upload registry

        run_id = uuid.uuid4().hex[:12]
        run_dir = RUNS_DIR / run_id
        run_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            "nextflow", "run", "main.nf",
            "-profile", req.profile,
            "-name", f"ws1-{run_id}",
            "-work-dir", str(run_dir / "work"),
            "--outdir", str(run_dir / "results"),
            "-with-weblog", f"{self.base_url}/runs/{run_id}/_weblog",
        ]
        if input_path:
            cmd += ["--input", input_path]
        if annotators:
            cmd += ["--annotators", ",".join(annotators)]

        info = RunInfo(
            id=run_id, status=RunStatus.running, profile=req.profile,
            annotators=annotators, input=input_path,
            started_at=time.time(), progress=TaskProgress(),
        )
        run = Run(info)
        self.runs[run_id] = run

        log = (run_dir / "nextflow.log").open("wb")
        proc = await asyncio.create_subprocess_exec(
            *cmd, cwd=str(PROJECT_DIR),
            stdout=log, stderr=asyncio.subprocess.STDOUT,
        )
        run.proc = proc
        asyncio.create_task(self._watch(run, proc, log))
        return info

    async def _watch(self, run: Run, proc, log) -> None:
        code = await proc.wait()
        log.close()
        run.info.exit_code = code
        if run.info.status != RunStatus.cancelled:
            run.info.status = RunStatus.completed if code == 0 else RunStatus.failed
        run.info.finished_at = time.time()
        run.publish({"event": "terminated", "status": run.info.status.value,
                     "exit_code": code})

    async def ingest_weblog(self, run_id: str, payload: dict) -> None:
        run = self.runs.get(run_id)
        if run is None:
            return
        event = payload.get("event")
        if event == "process_submitted":
            run.info.progress.submitted += 1
        elif event == "process_completed":
            run.info.progress.completed += 1
            if str((payload.get("trace") or {}).get("status")) == "FAILED":
                run.info.progress.failed += 1
        run.publish(payload)

    async def cancel(self, run_id: str) -> RunInfo | None:
        run = self.runs.get(run_id)
        if run is None:
            return None
        if run.proc is not None and run.proc.returncode is None:
            run.info.status = RunStatus.cancelled
            run.proc.terminate()
        return run.info

    # ---- reads ----
    def get(self, run_id: str) -> RunInfo | None:
        run = self.runs.get(run_id)
        return run.info if run else None

    def list(self) -> list[RunInfo]:
        return [r.info for r in self.runs.values()]

    def results(self, run_id: str) -> list[str]:
        # run_id is only a dict key; the path comes from the trusted run.dir.
        run = self.runs.get(run_id)
        if run is None:
            return []
        base = run.dir / "results"
        if not base.exists():
            return []
        return sorted(str(p.relative_to(base)) for p in base.rglob("*") if p.is_file())

    def result_path(self, run_id: str, rel: str) -> Path | None:
        run = self.runs.get(run_id)
        if run is None:
            return None
        listing = self.results(run_id)              # trusted relative paths (from rglob)
        if rel not in listing:
            return None                             # only files in the run's own listing
        # Build from the trusted listing entry, not the raw request string.
        return run.dir / "results" / listing[listing.index(rel)]
