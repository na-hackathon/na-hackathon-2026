"""Tests for bin/ws1-nftest (nf-test test with its log in .nf-test/).

nf-test is replaced by a script that records its arguments.
Run with:  python -m pytest tests/test_ws1_nftest.py
"""
from __future__ import annotations

import os
import subprocess
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin"


def test_nftest_wrapper_puts_the_log_in_the_work_dir(tmp_path):
    bindir, log = tmp_path / "fakebin", tmp_path / "calls.log"
    bindir.mkdir()
    (bindir / "nf-test").write_text(f'#!/usr/bin/env bash\necho "nf-test $*" >> {log}\n')
    (bindir / "nf-test").chmod(0o755)
    env = {**os.environ, "PATH": f"{bindir}:{os.environ['PATH']}"}
    subprocess.run([str(BIN / "ws1-nftest"), "tests/parse.nf.test", "--tag", "conda"], cwd=tmp_path, env=env, check=True)
    assert "nf-test test --log .nf-test/nf-test.log tests/parse.nf.test --tag conda" in log.read_text()
    assert (tmp_path / ".nf-test").is_dir()
