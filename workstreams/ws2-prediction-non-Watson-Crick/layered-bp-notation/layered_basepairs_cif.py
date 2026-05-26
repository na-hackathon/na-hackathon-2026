"""Layered base-pair notation built from a single mmCIF.

    _ndb_base_pair_list        - who pairs with whom (auth numbering + operator)
    _ndb_base_pair_annotation  - the Leontis-Westhof family as a string (l-w_family)
    _pdbx_struct_oper_list     - maps the operator id to its 1_555-style name

The two pair categories are joined on base_pair_id. The sequence still comes
from pdbx_poly_seq_scheme (or _atom_site as a fallback). Each base is keyed by
(chain, symmetry, number) so symmetry copies render as separate strands.

Usage:
    python3 layered_basepairs_cif.py <cif> [chains...] [--name NAME] [--block N]
                                     [--compact]

With no chains, all nucleic-acid chains that form pairs are used automatically,
so the CIF is the only required input.

Examples:
    python3 layered_basepairs_cif.py 9cfn.cif --name 9CFN     # auto chains
    python3 layered_basepairs_cif.py 9cfn.cif A --name 9CFN   # one chain

The notation prints to stdout; the round-trip check (True = lossless) to stderr.
"""

import argparse
import shlex
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


def _read_oper_map(lines: list[str]) -> dict[str, str]:
    """Map a struct_oper id to its symmetry name (e.g. '1' -> '1_555') from
    _pdbx_struct_oper_list, which appears either as a single key-value block or
    as a loop_."""
    kv = [l for l in lines if l.strip().startswith("_pdbx_struct_oper_list.")]
    if kv and all(len(l.split()) >= 2 for l in kv):       # single-row key-value form
        d = {l.split()[0].split(".")[1]: l.split()[1] for l in kv}
        return {d.get("id", "1"): d.get("name", IDENTITY)}
    idx, i = _loop_columns(lines, "_pdbx_struct_oper_list.")
    m: dict[str, str] = {}
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = shlex.split(row)              # quote-aware: 'crystal symmetry operation'
        if len(f) >= len(idx):
            m[f[idx["id"]]] = f[idx["name"]]
        i += 1
    return m


def _read_annotation(lines: list[str]) -> dict[str, str]:
    """Map base_pair_id -> Leontis-Westhof family string (cWW, tWH, ...) from
    _ndb_base_pair_annotation."""
    idx, i = _loop_columns(lines, "_ndb_base_pair_annotation.")
    c_id, c_lw = idx["base_pair_id"], idx["l-w_family"]
    lw_by_id: dict[str, str] = {}
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = row.split()
        if len(f) >= len(idx):
            lw_by_id[f[c_id]] = f[c_lw]
        i += 1
    return lw_by_id


def read_pairs(cif: Path, chains: set[str]) -> list[tuple[tuple, tuple, str]]:
    """De-duplicated pairs (Ra, Rb, lw) read from the CIF's own annotation, where
    each residue R is (chain, symmetry, number). _ndb_base_pair_list gives the
    residues (auth numbering) and operator ids; _ndb_base_pair_annotation gives
    the LW family, joined on base_pair_id. Ordering is canonical so the LW
    direction stays consistent."""
    lines = cif.read_text().splitlines()
    oper = _read_oper_map(lines)
    lw_by_id = _read_annotation(lines)
    idx, i = _loop_columns(lines, "_ndb_base_pair_list.")
    seen: dict[frozenset, tuple] = {}
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = row.split()
        if len(f) >= len(idx):
            ch1, ch2 = f[idx["auth_asym_id_1"]], f[idx["auth_asym_id_2"]]
            if ch1 in chains and ch2 in chains:
                s1 = oper.get(f[idx["struct_oper_id_1"]], IDENTITY)
                s2 = oper.get(f[idx["struct_oper_id_2"]], IDENTITY)
                Ra = (ch1, s1, int(f[idx["auth_seq_id_1"]]))
                Rb = (ch2, s2, int(f[idx["auth_seq_id_2"]]))
                lw = lw_by_id.get(f[idx["base_pair_id"]])
                if Ra > Rb:                # canonical order; flip the label to match
                    Ra, Rb, lw = Rb, Ra, _flip_lw(lw)
                seen[frozenset((Ra, Rb))] = (Ra, Rb, lw)
        i += 1
    return sorted(seen.values())


def list_chains(cif: Path) -> list[str]:
    """Chains that appear in the base-pair list, sorted. Used as the default when
    no chains are given on the command line; these are exactly the nucleic-acid
    chains that form annotated pairs, so protein chains are never included."""
    lines = cif.read_text().splitlines()
    idx, i = _loop_columns(lines, "_ndb_base_pair_list.")
    chains: set[str] = set()
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = row.split()
        if len(f) >= len(idx):
            chains.add(f[idx["auth_asym_id_1"]])
            chains.add(f[idx["auth_asym_id_2"]])
        i += 1
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


def _res_id(R: tuple) -> str:
    """Compact residue label: chain+number, with the operator only when it is a
    symmetry mate (e.g. 'A24', or 'A1012[2_755]')."""
    ch, sym, num = R
    return f"{ch}{num}" if sym == IDENTITY else f"{ch}{num}[{sym}]"


def _compact_line(pairs: list) -> str:
    """A sparse layer as an explicit pair list, e.g. 'A24,A31 A25,A29'."""
    return " ".join(f"{_res_id(Ra)},{_res_id(Rb)}" for Ra, Rb, _ in sorted(pairs))


def build(cif: Path, chains: list[str], name: str = "RNA",
          block: int | None = None, compact: bool = False) -> tuple[str, bool]:
    """Build the layered notation. Returns the text and the round-trip result.

    Args:
        cif: mmCIF file; provides both the sequence/numbering and the base pairs.
        chains: chain ids to render, in the order to lay them on the ruler.
        name: label for the header line.
        block: if set, wrap the lines into blocks of this many columns.
        compact: print the sparse non-cWW layers as explicit pair lists
            (e.g. 'A24,A31') instead of full-width dot-bracket; cWW stays
            dot-bracket. Saves space on large RNA, where most non-WC lines are
            almost all dots.
    """
    if not chains:
        chains = list_chains(cif)
    per_chain = read_residues(cif, chains)
    pairs = read_pairs(cif, set(chains))

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

    if compact:        # cWW stays dot-bracket; sparse non-cWW layers become pair lists
        rows = [f"{'seq':12}: " + seq_line]
        for (label, layer_pairs), (_, dbstr) in zip(layers, rendered):
            line = dbstr if label.split()[1] == "cWW" else _compact_line(layer_pairs)
            rows.append(f"{label:12}: " + line)
        text = head + "\n" + "\n".join(rows)
    elif block:
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
        description="Layered base-pair notation from a single annotated mmCIF "
                    "(no external pair list).")
    ap.add_argument("cif", type=Path,
                    help="mmCIF file with _ndb_base_pair_list/_annotation categories")
    ap.add_argument("chains", nargs="*",
                    help="chain ids, in the order to lay them on the ruler; "
                         "default: all nucleic-acid chains that form pairs")
    ap.add_argument("--name", default="RNA", help="label for the header line")
    ap.add_argument("--block", type=int, nargs="?", const=BLOCK, default=None,
                    help=f"wrap lines into blocks of this many columns "
                         f"(default {BLOCK} when given with no value)")
    ap.add_argument("--compact", action="store_true",
                    help="print sparse non-cWW layers as explicit pair lists "
                         "(e.g. A24,A31) instead of dot-bracket; saves space on "
                         "large RNA")
    a = ap.parse_args()
    text, ok = build(a.cif, a.chains, name=a.name, block=a.block, compact=a.compact)
    print(text)
    print(f"\n# round-trip recovers all pairs exactly: {ok}", file=sys.stderr)
    sys.exit(0 if ok else 1)
