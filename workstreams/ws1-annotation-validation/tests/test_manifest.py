"""Tests for bin/ws1-manifest (the run manifest).

Run with:  python -m pytest tests/test_manifest.py
"""
from __future__ import annotations

import hashlib
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin"
FIXTURES = Path(__file__).resolve().parents[1] / "tests" / "fixtures"
RUN_INFO = {"nextflow": "26.10.4", "command_line": "nextflow run main.nf -profile test", "run_name": "x",
            "profile": "test", "start": "2026-10-28T10:00:00Z", "params": {"outdir": "results"}}


def project(tmp_path, init_git=True):
    """A pipeline directory with envs/ws1.yml, committed when init_git."""
    directory = tmp_path / "proj"
    (directory / "envs").mkdir(parents=True)
    (directory / "envs" / "ws1.yml").write_text("name: ws1\n")
    if init_git:
        git = ["git", "-C", str(directory), "-c", "user.name=t", "-c", "user.email=t@t"]
        subprocess.run([*git[:3], "init", "-q"], check=True)
        subprocess.run([*git, "add", "."], check=True)
        subprocess.run([*git, "commit", "-q", "-m", "x"], check=True)
    return directory


def run_manifest(tmp_path, directory):
    structure = tmp_path / "s.cif"
    structure.write_text("data_x\n")
    run_info = tmp_path / "run_info.json"
    run_info.write_text(json.dumps(RUN_INFO))
    out = tmp_path / "run_manifest.json"
    subprocess.run([sys.executable, str(BIN / "ws1-manifest"), "--input", str(structure), "--run-info", str(run_info),
                    "--project-dir", str(directory), "--out", str(out), "--provenance",
                    str(FIXTURES / "rnapolis_provenance.json"), str(FIXTURES / "fr3d_provenance.json")], check=True)
    return json.loads(out.read_text())


@pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")
def test_manifest_collects_pipeline_input_annotators_and_environment(tmp_path):
    directory = project(tmp_path)
    head = subprocess.run(["git", "-C", str(directory), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    data = run_manifest(tmp_path, directory)
    assert data["pipeline"] == {"commit": head, "modified": False, "env_file": "envs/ws1.yml",
                                "env_file_sha256": hashlib.sha256(b"name: ws1\n").hexdigest()}
    assert data["run"] == RUN_INFO
    assert data["input"] == {"name": "s.cif", "sha256": hashlib.sha256(b"data_x\n").hexdigest()}
    assert [a["tool"] for a in data["annotators"]] == ["fr3d", "rnapolis"]            # sorted by file name
    assert data["environment"]["python"] == platform.python_version()
    assert "pytest" in data["environment"]["python_distributions"]


@pytest.mark.skipif(shutil.which("git") is None, reason="git is not installed")
def test_manifest_flags_uncommitted_changes(tmp_path):
    directory = project(tmp_path)
    (directory / "envs" / "ws1.yml").write_text("name: changed\n")
    assert run_manifest(tmp_path, directory)["pipeline"]["modified"] is True


def test_manifest_without_a_git_checkout_has_null_commit(tmp_path):
    pipeline = run_manifest(tmp_path, project(tmp_path, init_git=False))["pipeline"]
    assert pipeline["commit"] is None and pipeline["modified"] is None
