"""Unit tests for bin/comparison_table.py.

Run with:  python -m pytest tests/test_comparison_table.py
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parents[1] / "bin"


def make_run(root: Path, name: str, skipped=(0, 2)):
    run = root / name
    (run / "validation").mkdir(parents=True)
    (run / "basepairs").mkdir()
    (run / "validation" / "validation_report.json").write_text(json.dumps({
        "tools": ["fr3d", "rnapolis"], "model": "1",
        "pair_counts": {"fr3d": 10, "rnapolis": 12}, "total_distinct_pairs": 14,
        "shared_by_all_tools": 8, "agreed_by_all_tools": 7, "discordant_family": 1,
        "only_in": {"fr3d": 2, "rnapolis": 4}, "comp_id_mismatch": 0}))
    for tool, n in zip(("fr3d", "rnapolis"), skipped):
        (run / "basepairs" / f"{tool}_basepairs.skipped.tsv").write_text(
            "reason\tr1\tfamily\tr2\tdetail\n" + "x\ta\tcWW\tb\t\n" * n)
    return run


def run_script(*args):
    return subprocess.run([sys.executable, str(BIN / "comparison_table.py"), *map(str, args)],
                          capture_output=True, text=True)


def test_one_row_per_run_with_per_tool_columns(tmp_path):
    r = run_script(make_run(tmp_path, "8xzn"), make_run(tmp_path, "9cfn", skipped=(3, 0)), "--format", "tsv")
    assert r.returncode == 0, r.stderr
    header, first, second = [line.split("\t") for line in r.stdout.splitlines()]
    assert header == ["case", "model", "fr3d pairs", "rnapolis pairs", "fr3d skipped", "rnapolis skipped",
                      "distinct", "shared", "agreed", "discordant", "fr3d only", "rnapolis only", "name mismatch"]
    assert first == ["8xzn", "1", "10", "12", "0", "2", "14", "8", "7", "1", "2", "4", "0"]
    assert second[0] == "9cfn" and second[4:6] == ["3", "0"]


def test_markdown_is_the_default_format(tmp_path):
    lines = run_script(make_run(tmp_path, "8xzn")).stdout.splitlines()
    assert lines[0].startswith("| case | model |") and set(lines[1]) <= set("|-")


def test_single_annotator_report_is_rejected(tmp_path):
    run = make_run(tmp_path, "8xzn")
    (run / "validation" / "validation_report.json").write_text(json.dumps(
        {"tools": ["fr3d"], "model": "all", "pair_counts": {"fr3d": 5}}))
    r = run_script(run)
    assert r.returncode != 0 and "no cross-tool comparison" in r.stderr
