#!/usr/bin/env python3
"""Tabulate the validation reports of several pipeline runs, one row per run.

    comparison_table.py <outdir> [<outdir> ...] [--format md|tsv]

Each <outdir> is the --outdir of a run with more than one annotator; it provides
validation/validation_report.json and basepairs/<tool>_basepairs.skipped.tsv. A row is named
after its directory. Columns: model compared; per tool the distinct pairs in that model, the
rows skipped by the exporter and the pairs found by that tool only; then distinct pairs over all
tools, pairs shared by all tools, pairs with the same family in all tools (agreed), pairs shared
with different families (discordant) and shared pairs with different residue names.
"""
import argparse
import json
import sys
from pathlib import Path


def skipped_rows(outdir: Path, tool: str):
    path = outdir / "basepairs" / f"{tool}_basepairs.skipped.tsv"
    if not path.exists():
        return ""
    return str(len(path.read_text().splitlines()) - 1)       # minus the header


def table(outdirs):
    reports = {d: json.loads((d / "validation" / "validation_report.json").read_text()) for d in outdirs}
    tools = sorted({t for r in reports.values() for t in r["tools"]})
    header = (["case", "model"] + [f"{t} pairs" for t in tools] + [f"{t} skipped" for t in tools]
              + ["distinct", "shared", "agreed", "discordant"] + [f"{t} only" for t in tools]
              + ["name mismatch"])
    rows = []
    for d, r in reports.items():
        if "total_distinct_pairs" not in r:
            sys.exit(f"[comparison_table] {d}: no cross-tool comparison in the report (single annotator?)")
        rows.append([d.name, r["model"]]
                    + [str(r["pair_counts"].get(t, "")) for t in tools]
                    + [skipped_rows(d, t) for t in tools]
                    + [str(r["total_distinct_pairs"]), str(r["shared_by_all_tools"]),
                       str(r["agreed_by_all_tools"]), str(r["discordant_family"])]
                    + [str(r["only_in"].get(t, "")) for t in tools]
                    + [str(r["comp_id_mismatch"])])
    return header, rows


def render(header, rows, fmt):
    if fmt == "tsv":
        return "\n".join("\t".join(r) for r in [header] + rows)
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="Tabulate validation reports of several pipeline runs")
    ap.add_argument("outdirs", nargs="+", type=Path, help="--outdir of each run")
    ap.add_argument("--format", choices=["md", "tsv"], default="md")
    args = ap.parse_args()
    print(render(*table(args.outdirs), args.format))


if __name__ == "__main__":
    main()
