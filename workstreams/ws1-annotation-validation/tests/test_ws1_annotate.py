"""Tests for the checks of bin/ws1-annotate.

Run with:  python -m pytest tests/test_ws1_annotate.py
"""
from __future__ import annotations

import subprocess
from pathlib import Path

ANNOTATE = Path(__file__).resolve().parents[1] / "bin" / "ws1-annotate"


def annotate(tmp_path, *args):
    return subprocess.run(["bash", str(ANNOTATE), *args], capture_output=True, text=True)


def test_missing_submodule_is_named_with_the_fix(tmp_path):
    result = annotate(tmp_path, "--tool", "fr3d", "--structure", "s.cif", "--out", "o", "--src", str(tmp_path))
    assert result.returncode == 2
    assert "NA_pairwise_interactions.py not found" in result.stderr and "git submodule update --init" in result.stderr


def test_unsupported_tool_is_rejected(tmp_path):
    result = annotate(tmp_path, "--tool", "nope", "--structure", "s.cif", "--out", "o", "--src", str(tmp_path))
    assert result.returncode == 2 and "no annotator wired for tool 'nope'" in result.stderr
