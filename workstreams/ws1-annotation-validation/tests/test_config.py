"""Resolved nextflow.config: conda cache location and env file.

Run with:  python -m pytest tests/test_config.py   (skipped without nextflow on PATH)
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

WS1 = Path(__file__).resolve().parents[1]
pytestmark = pytest.mark.skipif(shutil.which("nextflow") is None, reason="nextflow is not installed")


def resolved(**env):
    environment = {k: v for k, v in os.environ.items() if k != "NXF_CONDA_CACHEDIR"}
    environment.update(env, NXF_ANSI_LOG="false", NXF_OFFLINE="true")
    out = subprocess.run(["nextflow", "config", "-flat"], cwd=WS1, env=environment, capture_output=True, text=True, check=True)
    return dict(line.split(" = ", 1) for line in out.stdout.splitlines() if " = " in line)


def test_conda_cache_defaults_to_a_directory_outside_work(tmp_path):
    config = resolved(HOME=str(tmp_path))
    assert config["conda.cacheDir"] == f"'{tmp_path}/.cache/ws1-conda'"


def test_conda_cache_follows_nxf_conda_cachedir(tmp_path):
    assert resolved(NXF_CONDA_CACHEDIR=str(tmp_path / "c"))["conda.cacheDir"] == f"'{tmp_path}/c'"


def test_env_file_defaults_to_the_single_env():
    assert resolved()["params.env_file"].strip("'") == str(WS1 / "envs" / "ws1.yml")
