"""Layered base-pair notation for a multi-chain complex, symmetry-aware.

Chains are laid end to end on one ruler, separated by '&' (the ViennaRNA cofold
convention), and a pair between two chains is just a bracket that opens in one
chain and closes in another. Each base is keyed by (chain, symmetry, number):
the symmetry field distinguishes a residue from a symmetry-generated copy of the
same chain, so self-complementary duplexes and crystal-packing contacts (where
FR3D reports a pair as chain A to chain A with a symmetry operator on one
partner) are represented as two strands rather than a spurious self-pairing.

The sequence and residue numbers come from pdbx_poly_seq_scheme when present,
otherwise from _atom_site (the coordinates), so structures lacking that category
(e.g. modeling results) still work.

Usage:
    python3 layered_basepairs.py <cif> <tsv> [chains...] [--name NAME] [--block N]

With no chains, all chains present in the FR3D TSV are used automatically.

Examples (mmCIF from files.rcsb.org/download/<ID>.cif, basepairs TSV from FR3D):
    # single chain, non-WC pairs and multi-pairing
    python3 layered_basepairs.py 9CFN.cif 9cfn_fr3d_basepairs.tsv A --name 9CFN
    # no chains given -> all chains in the TSV are used
    python3 layered_basepairs.py 9CFN.cif 9cfn_fr3d_basepairs.tsv --name 9CFN
    # pairs across chains (1XPE kissing-loop dimer)
    python3 layered_basepairs.py 1XPE.cif 1xpe_fr3d_basepairs.tsv A B --name 1XPE
    # crystal-symmetry duplex; the symmetry copy is picked up automatically
    python3 layered_basepairs.py 2Q1R.cif 2q1r_fr3d_basepairs.tsv A --name 2Q1R

The notation prints to stdout; the round-trip check (True = lossless) to stderr.
"""

import argparse
import sys
from pathlib import Path

# Bracket classes for crossing pairs within one layer (pseudoknot levels).
BRACKETS = ["()", "[]", "{}", "<>", "Aa", "Bb", "Cc", "Dd"]
OPENERS = set("([{<ABCD")
CLOSE_TO_OPEN = {")": "(", "]": "[", "}": "{", ">": "<",
                 "a": "A", "b": "B", "c": "C", "d": "D"}


def _flip_lw(lw: str) -> str:
    """Flip an LW code direction (cWH <-> cHW) when the pair is reordered."""
    if lw and len(lw) == 3 and lw[0] in "ct":
        return lw[0] + lw[2] + lw[1]
    return lw


def _cross(p: tuple[int, int], q: tuple[int, int]) -> bool:
    (i, j), (k, l) = p, q
    return i < k < j < l or k < i < l < j


# Fixed layer order: the 18 directed LW types, cWW first.
DIRECTED = ["cWW", "cWH", "cWS", "cHW", "cHH", "cHS", "cSW", "cSH", "cSS",
            "tWW", "tWH", "tWS", "tHW", "tHH", "tHS", "tSW", "tSH", "tSS"]

SEP = "&"          # chain boundary marker
BLOCK = 100        # default columns per block when --block is given with no value
IDENTITY = "1_555" # crystallographic identity operator (the asymmetric-unit copy)

# Map residue names to a single sequence letter (RNA and DNA); modified
# residues fall back to their last character, lowercased, so columns stay aligned.
NT = {"A": "A", "C": "C", "G": "G", "U": "U", "T": "T",
      "DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U"}


def _loop_columns(lines: list[str], prefix: str, start: int = 0) -> tuple[dict, int]:
    """Find a `prefix` loop, return ({column_name: index}, index of first data row)."""
    i = start
    while i < len(lines) and not lines[i].startswith(prefix):
        i += 1
    cols: list[str] = []
    while i < len(lines) and lines[i].startswith(prefix):
        cols.append(lines[i].strip().split(".")[1])
        i += 1
    return {n: k for k, n in enumerate(cols)}, i


def _read_poly_seq_scheme(lines: list[str], chains: list[str]) -> dict[str, list]:
    """Per-chain [(number, letter), ...] from pdbx_poly_seq_scheme."""
    idx, i = _loop_columns(lines, "_pdbx_poly_seq_scheme.")
    c_mon, c_num, c_strand = idx["mon_id"], idx["pdb_seq_num"], idx["pdb_strand_id"]
    per_chain: dict[str, list] = {c: [] for c in chains}
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "_", "loop_")):
            break
        f = row.split()
        if len(f) >= len(idx) and f[c_strand] in per_chain:
            mon = f[c_mon]
            per_chain[f[c_strand]].append((int(f[c_num]), NT.get(mon, mon[-1].lower())))
        i += 1
    return per_chain


def _read_atom_site(lines: list[str], chains: list[str]) -> dict[str, list]:
    """Per-chain [(number, letter), ...] from _atom_site (model 1, polymer only)."""
    idx, i = _loop_columns(lines, "_atom_site.")
    c_ch, c_num, c_comp = idx["auth_asym_id"], idx["auth_seq_id"], idx["label_comp_id"]
    c_lab = idx["label_seq_id"]
    c_model = idx.get("pdbx_PDB_model_num")
    per_chain: dict[str, list] = {c: [] for c in chains}
    seen: set[tuple[str, str]] = set()
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = row.split()
        if len(f) >= len(idx):
            if (c_model is not None and f[c_model] != "1") or f[c_lab] == ".":
                i += 1
                continue
            ch = f[c_ch]
            key = (ch, f[c_num])
            if ch in per_chain and key not in seen:
                seen.add(key)
                comp = f[c_comp]
                per_chain[ch].append((int(f[c_num]), NT.get(comp, comp[-1].lower())))
        i += 1
    return per_chain


def read_residues(cif: Path, chains: list[str]) -> dict[str, list]:
    """Per-chain [(number, letter), ...] from pdbx_poly_seq_scheme, or _atom_site
    as a fallback when that category is absent."""
    lines = cif.read_text().splitlines()
    has_scheme = any(l.startswith("_pdbx_poly_seq_scheme.") for l in lines)
    return (_read_poly_seq_scheme(lines, chains) if has_scheme
            else _read_atom_site(lines, chains))


def _sym(unit: list[str]) -> str:
    """Symmetry operator from an FR3D unit id. A missing operator denotes the
    asymmetric-unit copy, which is the crystallographic identity, so it is
    reported explicitly as 1_555. FR3D unit id is
    pdb|model|chain|comp|num|atom|alt|ins|symmetry."""
    return unit[8] if len(unit) > 8 and unit[8] else IDENTITY


def read_pairs(tsv: Path, chains: set[str]) -> list[tuple[tuple, tuple, str]]:
    """De-duplicated pairs (Ra, Rb, lw), where each residue R is (chain, symmetry,
    number). Both ends must be among `chains`; ordering is canonical so the LW
    direction is consistent."""
    seen: dict[frozenset, tuple] = {}
    for line in tsv.read_text().splitlines():
        c = line.strip().split("\t")
        if len(c) < 3:
            continue
        a, b = c[0].split("|"), c[2].split("|")
        if a[2] not in chains or b[2] not in chains:
            continue
        Ra, Rb, lw = (a[2], _sym(a), int(a[4])), (b[2], _sym(b), int(b[4])), c[1]
        if Ra > Rb:                       # canonical order; flip the label to match
            Ra, Rb, lw = Rb, Ra, _flip_lw(lw)
        seen[frozenset((Ra, Rb))] = (Ra, Rb, lw)
    return sorted(seen.values())


def list_chains(tsv: Path) -> list[str]:
    """Chains that appear in the FR3D pair list, sorted. Used as the default when
    no chains are given on the command line."""
    chains: set[str] = set()
    for line in tsv.read_text().splitlines():
        c = line.strip().split("\t")
        if len(c) < 3:
            continue
        a, b = c[0].split("|"), c[2].split("|")
        if len(a) > 2 and len(b) > 2:
            chains.add(a[2])
            chains.add(b[2])
    return sorted(chains)


def _pack(pairs: list) -> list[list]:
    """Pack same-type overflow pairs into as few layers as possible, each residue
    used at most once per layer."""
    buckets: list[list] = []
    for p in sorted(pairs):
        Ra, Rb, _ = p
        for bucket in buckets:
            keys = {k for x, y, _ in bucket for k in (x, y)}
            if Ra not in keys and Rb not in keys:
                bucket.append(p)
                break
        else:
            buckets.append([p])
    return buckets


def render(layer: list, key_to_col: dict, width: int, sep_cols: set[int]) -> str:
    """Render one layer to a dot-bracket line over the multi-strand ruler."""
    chars = [SEP if i in sep_cols else "." for i in range(width)]
    colpairs = [tuple(sorted((key_to_col[Ra], key_to_col[Rb]))) for Ra, Rb, _ in layer]
    classed: list[tuple[tuple[int, int], int]] = []
    for (i, j) in sorted(colpairs):
        cls = 0
        while any(c == cls and _cross((i, j), pq) for pq, c in classed):
            cls += 1
        classed.append(((i, j), cls))
        chars[i], chars[j] = BRACKETS[cls][0], BRACKETS[cls][1]
    return "".join(chars)


def parse(s: str, col_to_key: dict) -> set[frozenset]:
    """Read one rendered line back into {frozenset(Ra, Rb), ...}. '&' is skipped."""
    stacks: dict[str, list[int]] = {}
    out: set[frozenset] = set()
    for col, ch in enumerate(s):
        if ch in OPENERS:
            stacks.setdefault(ch, []).append(col)
        elif ch in CLOSE_TO_OPEN:
            o = stacks[CLOSE_TO_OPEN[ch]].pop()
            out.add(frozenset((col_to_key[o], col_to_key[col])))
    return out


def _label(strand: tuple[str, str]) -> str:
    ch, sym = strand
    return f"{ch}[{sym}]"


def _span(residues: list, strand: tuple[str, str]) -> str:
    nums = [R[2] for R in residues if (R[0], R[1]) == strand]
    return f"{nums[0]}-{nums[-1]}" if nums else "absent"


def build(cif: Path, tsv: Path, chains: list[str], name: str = "RNA",
          block: int | None = None) -> tuple[str, bool]:
    """Build the layered notation. Returns the text and the round-trip result.

    Args:
        cif: mmCIF file, used for the sequence and residue numbers.
        tsv: FR3D base-pair list.
        chains: chain ids to render, in the order to lay them on the ruler;
            if empty, all chains in the TSV are used.
        name: label for the header line.
        block: if set, wrap the lines into blocks of this many columns.
    """
    if not chains:
        chains = list_chains(tsv)
    per_chain = read_residues(cif, chains)
    pairs = read_pairs(tsv, set(chains))

    # One strand per (chain, symmetry): the identity copy (1_555) of each requested
    # chain, plus any symmetry copies that appear in the pairs (deterministic order).
    sym_seen = sorted({R[:2] for Ra, Rb, _ in pairs for R in (Ra, Rb) if R[1] != IDENTITY})
    strands: list[tuple[str, str]] = []
    for c in chains:
        strands.append((c, IDENTITY))
        strands.extend(s for s in sym_seen if s[0] == c)

    columns: list[tuple] = []          # ("res", R, letter) | ("sep",)
    key_to_col: dict[tuple, int] = {}
    residues: list[tuple] = []
    for n, (ch, sym) in enumerate(strands):
        if n > 0:
            columns.append(("sep",))
        for num, letter in per_chain.get(ch, []):
            R = (ch, sym, num)
            key_to_col[R] = len(columns)
            columns.append(("res", R, letter))
            residues.append(R)
    width = len(columns)
    sep_cols = {i for i, col in enumerate(columns) if col[0] == "sep"}
    col_to_key = {i: col[1] for i, col in enumerate(columns) if col[0] == "res"}

    slot = {t: i + 1 for i, t in enumerate(DIRECTED)}
    by_type: dict[str, list] = {}
    for Ra, Rb, lw in pairs:
        by_type.setdefault(lw, []).append((Ra, Rb, lw))

    layers, overflow = [], []           # overflow: (type, bucket), packed per type
    for t in DIRECTED:
        if t not in by_type:
            continue
        placed, used, extra = [], set(), []
        for Ra, Rb, lw in sorted(by_type[t]):
            if Ra in used or Rb in used:        # residue already used on this layer
                extra.append((Ra, Rb, lw))
            else:
                placed.append((Ra, Rb, lw))
                used |= {Ra, Rb}
        layers.append((f"L{slot[t]} {t}", placed))
        overflow += [(t, bucket) for bucket in _pack(extra)]
    for k, (t, bucket) in enumerate(overflow):
        layers.append((f"L{19 + k} {t}", bucket))   # type kept in the label

    seq_line = "".join(col[2] if col[0] == "res" else SEP for col in columns)
    rendered = [(label, render(p, key_to_col, width, sep_cols)) for label, p in layers]

    spans = ", ".join(f"{_label(s)}({_span(residues, s)})" for s in strands)
    head = f">{name} complex: chains {spans}   ('{SEP}' separates chains)"

    if block:
        body = "\n\n".join(_wrap(seq_line, rendered, columns, block))
        text = head + "\n\n" + body
    else:
        rows = [f"{'seq':12}: " + seq_line]
        rows += [f"{label:12}: " + string for label, string in rendered]
        text = head + "\n" + "\n".join(rows)

    ok = _roundtrip(rendered, col_to_key, pairs)
    return text, ok


def _roundtrip(rendered: list, col_to_key: dict, pairs: list) -> bool:
    """Parse the written lines back to pairs and require the exact input set. The
    LW type is taken from each layer's label, not from the input pairs, so this
    also checks that the type is recoverable from the notation alone."""
    recovered = set()
    for label, s in rendered:
        lw = label.split()[1]            # the LW type named by the layer label
        for keyset in parse(s, col_to_key):
            Ra, Rb = sorted(keyset)
            recovered.add((Ra, Rb, lw))
    return recovered == set(pairs)


def _wrap(seq_line: str, rendered: list, columns: list, block: int) -> list[str]:
    """Split the seq + layer lines into fixed-width blocks of `block` columns."""
    blocks = []
    for s in range(0, len(seq_line), block):
        e = min(s + block, len(seq_line))
        res = [col for col in columns[s:e] if col[0] == "res"]
        if res:
            r0, r1 = res[0][1], res[-1][1]
            tag = f"{_label(r0[:2])}:{r0[2]} .. {_label(r1[:2])}:{r1[2]}"
        else:
            tag = "separator"
        rows = [f"# cols {s}-{e - 1}  ({tag})",
                f"{'seq':12}: " + seq_line[s:e]]
        rows += [f"{label:12}: " + string[s:e] for label, string in rendered]
        blocks.append("\n".join(rows))
    return blocks


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Layered base-pair notation for an RNA complex (symmetry-aware).")
    ap.add_argument("cif", type=Path, help="mmCIF file (sequence + residue numbers)")
    ap.add_argument("tsv", type=Path, help="FR3D basepairs TSV")
    ap.add_argument("chains", nargs="*",
                    help="chain ids, in the order to lay them on the ruler; "
                         "default: all chains present in the TSV")
    ap.add_argument("--name", default="RNA", help="label for the header line")
    ap.add_argument("--block", type=int, nargs="?", const=BLOCK, default=None,
                    help=f"wrap lines into blocks of this many columns "
                         f"(default {BLOCK} when given with no value)")
    a = ap.parse_args()
    text, ok = build(a.cif, a.tsv, a.chains, name=a.name, block=a.block)
    print(text)
    print(f"\n# round-trip recovers all pairs exactly: {ok}", file=sys.stderr)
    sys.exit(0 if ok else 1)
