#!/usr/bin/env python3
"""Per-tool base-pair stats and (optional) cross-tool comparison for WS1 validation.

    compare_basepairs.py --annotations <a.cif> [<b.cif> ...] --out <prefix> [--compare]
                         [--key label|auth] [--model first|all|N] [--asu-only]

Always writes per-tool pair-count stats. With --compare (which ws1-validate passes only
when more than one tool is present) it additionally matches base pairs across tools and
reports which pairs every tool agreed on (same residues AND same L-W family), which differ
in family (discordant) and which only one tool reports. The >1-tool gate lives in
ws1-validate, not here.

Pair identity (--key label, the default): the model plus, for each of the two residues,
label_asym_id, label_seq_id, insertion code and symmetry operator (order-independent; the
identity operator counts as "no operator"). Residue names are not part of the identity;
they are compared at matched pairs and disagreements are counted as comp_id_mismatch.
Alternate position is not part of it either: tools differ in whether they report it, so
rows that differ only by it share a key and are counted in duplicate_keys. One physical
contact can still be written with different symmetry operators, so the identity is not
canonical for symmetry-related pairs; --asu-only keeps only pairs inside the asymmetric
unit, for a like-for-like comparison with tools that cannot report symmetry mates.

Model scope (--model, comparison only): 'first' = the lowest model present in every tool's
output (so a tool that annotated only the first model is compared like-for-like with one that
annotated all), 'all', or one model number.

--key auth reproduces the earlier behaviour: identity = auth chain, auth number and residue
name (no model, insertion code or symmetry operator), rows sharing a key resolved by the
last row, and all models (--model then defaults to 'all').
"""
import argparse
import json
import re
import sys
from pathlib import Path

from parse_basepairs import identity_operator_id, parse_basepairs


def tool_of(path):
    m = re.match(r"(.+?)_basepairs\.cif$", Path(path).name)
    return m.group(1) if m else Path(path).stem


def residue_ref(p, side, key):
    if key == "auth":
        return (p[f"auth_asym_id_{side}"], p[f"auth_seq_id_{side}"], p[f"comp_id_{side}"])
    return (p[f"label_asym_id_{side}"], p[f"label_seq_id_{side}"],
            p[f"ins_code_{side}"], p[f"struct_oper_id_{side}"])


def pair_key(p, key="label"):
    """(order-independent key string, (name_lo, name_hi)) of one base pair."""
    (lo, lo_name), (hi, hi_name) = sorted(
        ((residue_ref(p, s, key), p[f"comp_id_{s}"]) for s in ("1", "2")), key=lambda r: r[0])
    if key == "auth":
        return f"{lo[0]}|{lo[1]}|{lo[2]}_{hi[0]}|{hi[1]}|{hi[2]}", (lo_name, hi_name)

    def text(r):
        return f"{r[0]}|{r[1]}{r[2]}" + (f"@{r[3]}" if r[3] else "")
    return f"m{p['model']}:{text(lo)}~{text(hi)}", (lo_name, hi_name)


def normalize_operators(pairs, identity):
    """Write the identity operator as '' so it compares equal across tools."""
    for p in pairs:
        for side in ("1", "2"):
            if p[f"struct_oper_id_{side}"] == identity:
                p[f"struct_oper_id_{side}"] = ""


def in_asu(p):
    return p["struct_oper_id_1"] == "" and p["struct_oper_id_2"] == ""


def choose_model(pairs_by_tool, spec):
    """Model number (str) to compare, or None for all models."""
    if spec == "all":
        return None
    if spec != "first":
        return spec
    present = [{p["model"] for p in ps} for ps in pairs_by_tool.values() if ps]
    if not present:
        return None
    common = set.intersection(*present)
    if not common:
        sys.exit("[compare_basepairs] no model is present in every tool's output; "
                 "use --model all or --model N")
    return min(common, key=lambda m: (not m.isdigit(), int(m) if m.isdigit() else 0, m))


def index_pairs(pairs, key):
    """key -> [(family, names), ...] in file order; rows sharing a key are kept, not overwritten."""
    idx = {}
    for p in pairs:
        k, names = pair_key(p, key)
        idx.setdefault(k, []).append((p["family"], names))
    return idx


def families(rows, key):
    fams = [f for f, _ in rows]
    return {fams[-1]} if key == "auth" else set(fams)   # --key auth keeps the earlier last-row-wins


def family_label(rows, key):
    return "/".join(sorted(families(rows, key))) if rows else "."


def summarize(per_tool, key, compare):
    """per_tool: {tool: index}. Returns (stats dict, TSV rows or None)."""
    tools = sorted(per_tool)
    rows = {t: sum(len(v) for v in per_tool[t].values()) for t in tools}
    stats = {"pair_counts": {t: len(per_tool[t]) for t in tools},
             "duplicate_keys": {t: rows[t] - len(per_tool[t]) for t in tools}}
    if not compare:
        return stats, None

    all_keys = sorted(set().union(*[set(d) for d in per_tool.values()]))
    shared = agreed = name_mismatch = 0
    only = {t: 0 for t in tools}
    table = []
    for k in all_keys:
        found = [t for t in tools if k in per_tool[t]]
        fams = {t: family_label(per_tool[t].get(k), key) for t in tools}
        same = False
        if len(found) == len(tools):
            shared += 1
            same = len({frozenset(families(per_tool[t][k], key)) for t in tools}) == 1
            agreed += 1 if same else 0
            if len({names for t in tools for _, names in per_tool[t][k]}) > 1:
                name_mismatch += 1
        elif len(found) == 1:
            only[found[0]] += 1
        table.append((k, fams, found, same))
    stats.update({
        "total_distinct_pairs": len(all_keys),
        "agreed_by_all_tools": agreed,      # present in every tool AND same L-W family
        "shared_by_all_tools": shared,
        "discordant_family": shared - agreed,
        "only_in": only,
        "comp_id_mismatch": name_mismatch,
    })
    return stats, table


def main():
    ap = argparse.ArgumentParser(description="WS1 base-pair stats + optional cross-tool comparison")
    ap.add_argument("--annotations", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--compare", action="store_true",
                    help="also run the cross-tool comparison (ws1-validate sets this when >1 tool)")
    ap.add_argument("--key", choices=["label", "auth"], default="label",
                    help="pair identity: label ids + insertion code + symmetry operator + model "
                         "(default), or the earlier auth chain/number/name key")
    ap.add_argument("--model", default=None,
                    help="comparison only: 'first' (default for --key label) = lowest model present in "
                         "every tool, 'all' (default for --key auth), or a model number")
    ap.add_argument("--asu-only", action="store_true",
                    help="keep only pairs whose two residues are both in the asymmetric unit")
    args = ap.parse_args()

    tool_files = {tool_of(p): p for p in args.annotations}
    if len(tool_files) < len(args.annotations):
        sys.exit("[compare_basepairs] two annotation files map to the same tool name")
    tools = sorted(tool_files)
    pairs = {t: parse_basepairs(tool_files[t]) for t in tools}
    for t in tools:
        normalize_operators(pairs[t], identity_operator_id(tool_files[t]))

    report = {"tools": tools, "key": args.key, "model": "all",
              "rows": {t: len(pairs[t]) for t in tools},
              "rows_with_symmetry_operator": {t: sum(1 for p in pairs[t] if not in_asu(p)) for t in tools}}

    if args.compare:
        model = choose_model(pairs, args.model or ("all" if args.key == "auth" else "first"))
        if model is not None:
            report["model"] = model
            report["rows_outside_compared_model"] = {
                t: sum(1 for p in pairs[t] if p["model"] != model) for t in tools}
            pairs = {t: [p for p in ps if p["model"] == model] for t, ps in pairs.items()}

    if args.asu_only:
        report["rows_dropped_non_asu"] = {t: sum(1 for p in pairs[t] if not in_asu(p)) for t in tools}
        pairs = {t: [p for p in ps if in_asu(p)] for t, ps in pairs.items()}

    stats, table = summarize({t: index_pairs(pairs[t], args.key) for t in tools}, args.key, args.compare)
    report.update(stats)

    if args.compare:
        header = ["pair"] + [f"family_{t}" for t in tools] + ["n_tools", "all_agree"]
        lines = ["\t".join(header)]
        for k, fams, found, same in table:
            lines.append("\t".join([k] + [fams[t] for t in tools]
                                   + [str(len(found)), "yes" if same else "no"]))
        msg = (f"{len(tools)} tools, {report['total_distinct_pairs']} distinct pairs, "
               f"{report['agreed_by_all_tools']} agreed by all "
               f"({report['discordant_family']} discordant, only_in={report['only_in']}, "
               f"comp_id_mismatch={report['comp_id_mismatch']})")
    else:
        report["comparison"] = "skipped (single tool)"
        lines = ["tool\tbase_pairs"] + [f"{t}\t{report['pair_counts'][t]}" for t in tools]
        msg = f"per-tool stats only ({', '.join(tools) or 'none'}); no cross-tool comparison"

    Path(f"{args.out}.json").write_text(json.dumps(report, indent=2) + "\n")
    Path(f"{args.out}.tsv").write_text("\n".join(lines) + "\n")
    print(f"[compare_basepairs] {msg}", file=sys.stderr)


if __name__ == "__main__":
    main()
