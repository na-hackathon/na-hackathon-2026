"""Tests for bin/ws1-env (conda env built once, interrupted builds redone, builds serialised).

conda is replaced by a script that creates the prefix and records the call.
Run with:  python -m pytest tests/test_ws1_env.py
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin"


def shim(directory: Path, name: str, body: str) -> None:
    script = directory / name
    script.write_text("#!/usr/bin/env bash\n" + body)
    script.chmod(0o755)


@pytest.fixture
def fake(tmp_path):
    """A PATH with a fake conda (creates the prefix, logs the call)."""
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    log = tmp_path / "calls.log"
    shim(bindir, "conda", f'echo "conda $*" >> {log}\n'
                          'prefix=; while [ $# -gt 0 ]; do [ "$1" = "-p" ] && prefix=$2; shift; done\n'
                          'sleep 0.3; mkdir -p "$prefix/bin"\n')
    env = {**os.environ, "PATH": f"{bindir}:{os.environ['PATH']}", "NXF_CONDA_CACHEDIR": str(tmp_path / "cache")}
    return env, log, tmp_path


def ws1_env(env, *args):
    return subprocess.run([str(BIN / "ws1-env"), *args], env=env, capture_output=True, text=True)


def conda_calls(log: Path) -> int:
    return len([line for line in log.read_text().splitlines() if line.startswith("conda ")]) if log.exists() else 0


def test_env_is_built_once_and_its_prefix_printed(fake):
    env, log, tmp = fake
    first, second = ws1_env(env), ws1_env(env)
    assert first.returncode == 0 and first.stdout == second.stdout
    prefix = Path(first.stdout.strip())
    assert prefix.parent == tmp / "cache" and prefix.name.startswith("ws1-") and (prefix / ".ws1-complete").exists()
    assert conda_calls(log) == 1                               # the second call found the finished env


def test_interrupted_build_is_removed_and_redone(fake):
    env, log, tmp = fake
    prefix = Path(ws1_env(env).stdout.strip())
    (prefix / ".ws1-complete").unlink()                        # what a killed build leaves behind
    (prefix / "half-built").write_text("x")
    assert ws1_env(env).returncode == 0
    assert (prefix / ".ws1-complete").exists() and not (prefix / "half-built").exists()
    assert conda_calls(log) == 2


def test_another_env_file_gets_another_prefix(fake):
    env, log, tmp = fake
    other = tmp / "other.yml"
    other.write_text("name: other\n")
    assert ws1_env(env).stdout != ws1_env(env, str(other)).stdout
    assert conda_calls(log) == 2


@pytest.mark.skipif(shutil.which("flock") is None, reason="flock is not installed")
def test_concurrent_calls_build_the_env_once(fake):
    env, log, tmp = fake
    processes = [subprocess.Popen([str(BIN / "ws1-env")], env=env, stdout=subprocess.PIPE, text=True) for _ in range(4)]
    prefixes = {p.communicate()[0].strip() for p in processes}
    assert len(prefixes) == 1 and all(p.returncode == 0 for p in processes)
    assert conda_calls(log) == 1


def test_failed_build_fails_and_leaves_no_completed_env(fake):
    env, log, tmp = fake
    shim(Path(env["PATH"].split(":")[0]), "conda", "exit 3\n")
    result = ws1_env(env)
    assert result.returncode != 0 and not list((tmp / "cache").glob("ws1-*/.ws1-complete"))

