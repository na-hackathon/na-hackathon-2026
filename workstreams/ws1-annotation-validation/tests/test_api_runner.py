"""Unit tests for the run-launching helpers of api/runner.py (no Nextflow process is started).

Needs pydantic (a dependency of the API); skipped where it is not installed.
Run with:  python -m pytest tests/test_api_runner.py
"""
from __future__ import annotations

import asyncio
import re
import sys
from pathlib import Path

import pytest

pytest.importorskip("pydantic")
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from api import runner  # noqa: E402
from api.models import RunRequest  # noqa: E402


def test_profiles_are_defined_in_nextflow_config():
    config = (runner.PROJECT_DIR / "nextflow.config").read_text()
    for profile in runner.PROFILES:
        assert re.search(rf"^\s+{profile}\s*\{{", config, re.M), profile


def test_annotators_are_wired_in_ws1_annotate():
    script = (runner.PROJECT_DIR / "bin" / "ws1-annotate").read_text()
    for annotator in runner.ANNOTATORS:
        assert f"\n    {annotator})" in script, annotator


def test_command_has_no_profile_by_default(tmp_path):
    cmd = runner.build_command("abc", tmp_path, "http://host:8000")
    assert cmd[:5] == ["nextflow", "-log", str(tmp_path / ".nextflow.log"), "run", "main.nf"]
    assert "-profile" not in cmd
    assert cmd[cmd.index("-work-dir") + 1] == str(tmp_path / "work")
    assert cmd[cmd.index("--outdir") + 1] == str(tmp_path / "results")
    assert cmd[cmd.index("-with-weblog") + 1] == "http://host:8000/runs/abc/_weblog"


def test_command_with_profile_input_and_annotators(tmp_path):
    cmd = runner.build_command("abc", tmp_path, "http://host:8000", "test", "/up/x.cif", ["fr3d", "rnapolis"])
    assert cmd[cmd.index("-profile") + 1] == "test"
    assert cmd[cmd.index("--input") + 1] == "/up/x.cif"
    assert cmd[cmd.index("--annotators") + 1] == "fr3d,rnapolis"


def test_command_passes_the_prebuilt_env(tmp_path):
    cmd = runner.build_command("abc", tmp_path, "http://host:8000", env_file="/cache/ws1-0123456789ab")
    assert cmd[cmd.index("--env_file") + 1] == "/cache/ws1-0123456789ab"
    assert "--env_file" not in runner.build_command("abc", tmp_path, "http://host:8000")


def script(path, body):
    path.write_text("#!/bin/sh\n" + body)
    path.chmod(0o755)
    return path


def test_ensure_env_returns_the_prefix_printed_by_ws1_env(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "UPLOAD_DIR", tmp_path / "uploads")
    monkeypatch.setattr(runner, "ENV_SCRIPT", script(tmp_path / "ws1-env", 'echo "/cache/ws1-abc"\n'))
    assert asyncio.run(runner.RunManager("http://host:8000").ensure_env()) == "/cache/ws1-abc"


def test_ensure_env_reports_a_failed_build(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "UPLOAD_DIR", tmp_path / "uploads")
    monkeypatch.setattr(runner, "ENV_SCRIPT", script(tmp_path / "ws1-env", 'echo "solver failed" >&2\nexit 1\n'))
    with pytest.raises(RuntimeError, match="solver failed"):
        asyncio.run(runner.RunManager("http://host:8000").ensure_env())


def test_launch_rejects_a_profile_missing_from_nextflow_config(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, "UPLOAD_DIR", tmp_path / "uploads")
    manager = runner.RunManager("http://host:8000")
    with pytest.raises(ValueError, match="unknown profile 'conda'"):
        asyncio.run(manager.launch(RunRequest(profile="conda")))
