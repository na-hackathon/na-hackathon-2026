"""Pydantic models for the WS1 control-plane API."""
from __future__ import annotations

from enum import Enum

from pydantic import BaseModel, Field


class RunStatus(str, Enum):
    running = "running"
    completed = "completed"
    failed = "failed"
    cancelled = "cancelled"


class UploadInfo(BaseModel):
    """Result of POST /uploads — an opaque id to pass back as a run's input."""
    upload_id: str
    filename: str


class RunRequest(BaseModel):
    """A request to launch a pipeline run."""
    input: str | None = Field(
        default=None,
        description="An upload_id from POST /uploads. Omit with profile=test, "
                    "which uses the bundled 9CFN.",
    )
    annotators: list[str] | None = Field(
        default=None, description="Annotator names; defaults to the pipeline's config."
    )
    profile: str = Field(default="conda", description="Nextflow -profile to use.")


class TaskProgress(BaseModel):
    submitted: int = 0
    completed: int = 0
    failed: int = 0


class RunInfo(BaseModel):
    id: str
    status: RunStatus
    profile: str
    annotators: list[str] | None = None
    input: str | None = None
    started_at: float
    finished_at: float | None = None
    exit_code: int | None = None
    progress: TaskProgress = TaskProgress()
