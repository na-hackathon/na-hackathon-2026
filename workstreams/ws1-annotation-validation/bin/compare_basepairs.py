#!/usr/bin/env python3
"""Per-tool base-pair stats and (optional) cross-tool comparison for WS1 validation.

    compare_basepairs.py --annotations <a.cif> [<b.cif> ...] --out <prefix> [--compare]

Always writes per-tool pair-count stats. With --compare (which ws1-validate passes only
when more than one tool is present) it additionally keys every base pair on its two
residues (order-independent) and reports which pairs every tool agreed on (same residues
AND same L-W family). The >1-tool gate lives in ws1-validate, not here.
"""
import argparse
import json
import re
import sys
from pathlib import Path

from parse_basepairs import parse_basepairs


def tool_of(path):
    m = re.match(r"(.+?)_basepairs\.cif$", Path(path).name)
    return m.group(1) if m else Path(path).stem


def pair_key(p):
    a = (p["auth_asym_id_1"], p["auth_seq_id_1"], p["comp_id_1"])
    b = (p["auth_asym_id_2"], p["auth_seq_id_2"], p["comp_id_2"])
    lo, hi = sorted([a, b])
    return f"{lo[0]}|{lo[1]}|{lo[2]}_{hi[0]}|{hi[1]}|{hi[2]}"


def main():
    ap = argparse.ArgumentParser(description="WS1 base-pair stats + optional cross-tool comparison")
    ap.add_argument("--annotations", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--compare", action="store_true",
                    help="also run the cross-tool comparison (ws1-validate sets this when >1 tool)")
    args = ap.parse_args()

    per_tool = {tool_of(p): {pair_key(x): x["family"] for x in parse_basepairs(p)}
                for p in args.annotations}
    tools = sorted(per_tool)
    report = {"tools": tools, "pair_counts": {t: len(per_tool[t]) for t in tools}}

    if args.compare:
        all_keys = sorted(set().union(*[set(d) for d in per_tool.values()]))
        rows, agreed = [], 0
        for k in all_keys:
            fams = {t: per_tool[t].get(k) for t in tools}
            found = [t for t in tools if fams[t] is not None]
            same = len(found) == len(tools) and len({fams[t] for t in found}) == 1
            agreed += 1 if same else 0
            rows.append((k, fams, found, same))
        report["total_distinct_pairs"] = len(all_keys)
        report["agreed_by_all_tools"] = agreed
        header = ["pair"] + [f"family_{t}" for t in tools] + ["n_tools", "all_agree"]
        lines = ["\t".join(header)]
        for k, fams, found, same in rows:
            lines.append("\t".join([k] + [fams[t] or "." for t in tools]
                                   + [str(len(found)), "yes" if same else "no"]))
        msg = f"{len(tools)} tools, {len(all_keys)} distinct pairs, {agreed} agreed by all"
    else:
        report["comparison"] = "skipped (single tool)"
        lines = ["tool\tbase_pairs"] + [f"{t}\t{report['pair_counts'][t]}" for t in tools]
        msg = f"per-tool stats only ({', '.join(tools) or 'none'}); no cross-tool comparison"

    Path(f"{args.out}.json").write_text(json.dumps(report, indent=2) + "\n")
    Path(f"{args.out}.tsv").write_text("\n".join(lines) + "\n")
    print(f"[compare_basepairs] {msg}", file=sys.stderr)


if __name__ == "__main__":
    main()
