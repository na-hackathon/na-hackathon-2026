"""Layered base-pair notation + 3D view, from a single DNATCO-extended mmCIF.

One command does everything:
  1. reads the base pairs straight from the CIF's own annotation,
  2. prints the layered dot-bracket notation (lossless; round-trip checked),
  3. prints iCn3D commands you paste into iCn3D to draw every pair on the real 3D
     structure, each pair its own colour (canonical Watson-Crick pairs grey).

You load the structure in iCn3D, paste the commands into its command box, press
Enter. No fragile long URLs. Use --script FILE to load them via File > Open instead.

Usage:
    python3 standalone_lbn_script.py <cif> [chains...] [--name NAME] [--id PDBID]
                                     [--noncanonical] [--compact] [--block N]
                                     [--script FILE]

    python3 standalone_lbn_script.py 9cfn.cif --name 9CFN          # auto chains
    python3 standalone_lbn_script.py 1XPE.cif A B --name 1XPE      # pick chains
    python3 standalone_lbn_script.py 9cfn.cif --noncanonical       # only non-canonical (notation + 3D)
    python3 standalone_lbn_script.py 9cfn.cif --script 9cfn.icn3d  # write a script file (any size)

The notation prints to stdout; the round-trip result and the iCn3D commands (or
a --script file path) follow.
"""

import argparse
import shlex
import sys
from pathlib import Path

# --------------------------------------------------------------------------- #
#  Constants                                                                   #
# --------------------------------------------------------------------------- #

# Bracket classes for crossing pairs within one layer (pseudoknot levels).
BRACKETS = ["()", "[]", "{}", "<>", "Aa", "Bb", "Cc", "Dd"]
OPENERS = set("([{<ABCD")
CLOSE_TO_OPEN = {")": "(", "]": "[", "}": "{", ">": "<",
                 "a": "A", "b": "B", "c": "C", "d": "D"}

# Fixed layer order: the 18 directed LW types, cWW first.
DIRECTED = ["cWW", "cWH", "cWS", "cHW", "cHH", "cHS", "cSW", "cSH", "cSS",
            "tWW", "tWH", "tWS", "tHW", "tHH", "tHS", "tSW", "tSH", "tSS"]

SEP = "&"           # chain boundary marker
BLOCK = 100         # default columns per block when --block is given with no value
IDENTITY = "1_555"  # crystallographic identity operator (asymmetric-unit copy)

# Residue name -> one sequence letter (RNA and DNA); modified residues fall back
# to their last character, lowercased, so columns stay aligned.
NT = {"A": "A", "C": "C", "G": "G", "U": "U", "T": "T",
      "DA": "A", "DC": "C", "DG": "G", "DT": "T", "DU": "U"}

ICN3D_URL = "https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html"

# Bold, well-separated colours so non-canonical pairs read from a distance.
# cWW (the canonical helix) is muted grey; every other family gets a vivid hue.
CWW_COLOR = "888888"
BOLD = ["FF00FF", "FF8000", "00C800", "00E5E5", "FFD000", "FF0000",
        "0000FF", "8000FF", "FF0080", "00FF80", "A0FF00", "0080FF",
        "FF4060", "8B4513", "FF80C0", "00B0FF", "9400D3"]


# --------------------------------------------------------------------------- #
#  Small helpers                                                               #
# --------------------------------------------------------------------------- #

def _flip_lw(lw: str) -> str:
    """Flip an LW code direction (cWH <-> cHW) when a pair is reordered."""
    if lw and len(lw) == 3 and lw[0] in "ct":
        return lw[0] + lw[2] + lw[1]
    return lw


def _cross(p: tuple[int, int], q: tuple[int, int]) -> bool:
    (i, j), (k, l) = p, q
    return i < k < j < l or k < i < l < j


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


# --------------------------------------------------------------------------- #
#  Read the sequence / numbering from the CIF                                  #
# --------------------------------------------------------------------------- #

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


# --------------------------------------------------------------------------- #
#  Read the base pairs from the CIF's own annotation                           #
# --------------------------------------------------------------------------- #

def _read_oper_map(lines: list[str]) -> dict[str, str]:
    """Map a struct_oper id to its symmetry name (e.g. '1' -> '1_555')."""
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
    """Map base_pair_id -> Leontis-Westhof family string from
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
    """De-duplicated pairs (Ra, Rb, lw) read from the CIF annotation, where each
    residue R is (chain, symmetry, number). _ndb_base_pair_list gives the residues
    and operators; _ndb_base_pair_annotation gives the LW family, joined on
    base_pair_id. Ordering is canonical so the LW direction stays consistent."""
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
    """Chains that appear in the base-pair list, sorted. Used when no chains are
    given; these are exactly the nucleic-acid chains that form annotated pairs."""
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


# --------------------------------------------------------------------------- #
#  Build the layered notation                                                  #
# --------------------------------------------------------------------------- #

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


def build_notation(cif: Path, chains: list[str], name: str = "RNA",
                   block: int | None = None, compact: bool = False,
                   noncanonical: bool = False) -> tuple[str, bool]:
    """Build the layered notation. Returns the text and the round-trip result.

    compact: print the sparse non-cWW layers as explicit pair lists (e.g.
    'A24,A31') instead of full-width dot-bracket; cWW stays dot-bracket.
    noncanonical: keep only non-canonical pairs (drop true Watson-Crick
    A-U/G-C/A-T); judged by the bases, so a cWW U-U/U-G wobble stays on L1."""
    if not chains:
        chains = list_chains(cif)
    per_chain = read_residues(cif, chains)
    pairs = read_pairs(cif, set(chains))

    if noncanonical:
        base_of = {(ch, num): letter
                   for ch, lst in per_chain.items() for num, letter in lst}
        pairs = [p for p in pairs
                 if not _is_canonical(p[2], base_of.get((p[0][0], p[0][2]), "?"),
                                            base_of.get((p[1][0], p[1][2]), "?"))]

    # One strand per (chain, symmetry): identity copy of each chain, plus any
    # symmetry copies that appear in the pairs (deterministic order).
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
    """Parse the written lines back to pairs and require the exact input set, with
    the LW type taken from each layer's label."""
    recovered = set()
    for label, s in rendered:
        lw = label.split()[1]
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


# --------------------------------------------------------------------------- #
#  Build the iCn3D 3D overlay (same pairs, plus coordinates from the CIF)      #
# --------------------------------------------------------------------------- #

WC_PAIRS = ({"A", "U"}, {"G", "C"}, {"A", "T"})   # canonical Watson-Crick base pairs


def _is_canonical(lw: str, base_a: str, base_b: str) -> bool:
    """A canonical Watson-Crick pair: cWW geometry AND A-U/G-C/A-T bases. cWW is a
    *family*, not 'canonical' -- a cWW U-U or G-U is non-canonical."""
    return lw == "cWW" and {base_a.upper(), base_b.upper()} in WC_PAIRS


def read_asu_c1(cif: Path, chains: set[str]) -> dict[tuple, tuple[float, float, float]]:
    """{(chain, num): (x, y, z)} for the asymmetric unit, using each residue's C1'
    atom (falls back to the residue's first atom if C1' is missing)."""
    lines = cif.read_text().splitlines()
    idx, i = _loop_columns(lines, "_atom_site.")
    c_ch, c_num, c_atom = idx["auth_asym_id"], idx["auth_seq_id"], idx["label_atom_id"]
    c_x, c_y, c_z = idx["Cartn_x"], idx["Cartn_y"], idx["Cartn_z"]
    c_model = idx.get("pdbx_PDB_model_num")
    coords: dict[tuple, tuple] = {}
    fallback: dict[tuple, tuple] = {}
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = row.split()
        if len(f) >= len(idx):
            if c_model is not None and f[c_model] != "1":
                i += 1
                continue
            ch = f[c_ch]
            if ch in chains:
                key = (ch, int(f[c_num]))
                xyz = (float(f[c_x]), float(f[c_y]), float(f[c_z]))
                if f[c_atom].strip('"') == "C1'":
                    coords[key] = xyz
                else:
                    fallback.setdefault(key, xyz)
        i += 1
    for k, v in fallback.items():
        coords.setdefault(k, v)
    return coords


def read_oper_matrices(cif: Path) -> dict[str, tuple[tuple, tuple]]:
    """{operator_name: (3x3 row-major R, 3-vector t)} from _pdbx_struct_oper_list.
    Matrices are already Cartesian, so a symmetry mate of point p is R*p + t."""
    lines = cif.read_text().splitlines()
    out: dict[str, tuple] = {IDENTITY: ((1., 0., 0., 0., 1., 0., 0., 0., 1.),
                                        (0., 0., 0.))}
    if not any(l.strip().startswith("_pdbx_struct_oper_list.matrix[1][1]") for l in lines):
        return out
    idx, i = _loop_columns(lines, "_pdbx_struct_oper_list.")
    rk = ("matrix[1][1]", "matrix[1][2]", "matrix[1][3]",
          "matrix[2][1]", "matrix[2][2]", "matrix[2][3]",
          "matrix[3][1]", "matrix[3][2]", "matrix[3][3]")
    while i < len(lines):
        row = lines[i].strip()
        if row.startswith(("#", "loop_", "_")):
            break
        f = shlex.split(row)
        if len(f) >= len(idx):
            out[f[idx["name"]]] = (tuple(float(f[idx[k]]) for k in rk),
                                   tuple(float(f[idx[k]]) for k in
                                         ("vector[1]", "vector[2]", "vector[3]")))
        i += 1
    return out


def _apply(R: tuple, asu: dict, opers: dict) -> tuple | None:
    """Cartesian C1' coordinate of residue R=(chain, symmetry, num), applying the
    symmetry operator to the asymmetric-unit position. None if unavailable."""
    chain, sym, num = R
    base = asu.get((chain, num))
    if base is None:
        return None
    if sym == IDENTITY:
        return base
    M = opers.get(sym)
    if M is None:
        return None
    r, t = M
    x, y, z = base
    return (r[0] * x + r[1] * y + r[2] * z + t[0],
            r[3] * x + r[4] * y + r[5] * z + t[1],
            r[6] * x + r[7] * y + r[8] * z + t[2])


def build_icn3d_lines(cif: Path, chains: list[str], noncanonical: bool = False
                      ) -> tuple[list[str], bool, int, dict[str, str]]:
    """iCn3D 'add line' commands for the same pairs the notation uses, drawn
    between C1' atoms. Every canonical Watson-Crick pair is the same thin grey;
    each non-canonical pair (including a cWW U-U) gets its own bold thick colour.
    Returns (commands, used a symmetry mate?, skipped count, legend list)."""
    pairs = read_pairs(cif, set(chains))
    asu = read_asu_c1(cif, set(chains))
    opers = read_oper_matrices(cif)
    base_of = {(ch, num): letter
               for ch, lst in read_residues(cif, list(chains)).items()
               for num, letter in lst}

    cmds: list[str] = []
    used_sym = False
    skipped = 0
    legend: list[tuple[str, str]] = []   # (colour, "pair family") per non-canonical pair
    nc = 0                               # running index -> a distinct colour per pair
    for Ra, Rb, lw in pairs:
        # Canonical = Watson-Crick base pair, judged by the bases (not the family),
        # so a cWW U-U is non-canonical and still gets drawn.
        ba = base_of.get((Ra[0], Ra[2]), "?")
        bb = base_of.get((Rb[0], Rb[2]), "?")
        canonical = _is_canonical(lw, ba, bb)
        if noncanonical and canonical:
            continue
        pa, pb = _apply(Ra, asu, opers), _apply(Rb, asu, opers)
        if pa is None or pb is None:
            skipped += 1
            continue
        if Ra[1] != IDENTITY or Rb[1] != IDENTITY:
            used_sym = True
        if canonical:
            color, radius = CWW_COLOR, "0.3"             # every canonical WC pair: grey, thin
        else:
            color, radius = BOLD[nc % len(BOLD)], "0.8"  # each non-canonical pair: own colour
            nc += 1
            legend.append((color, f"{lw} {ba}-{bb} {_res_id(Ra)}-{_res_id(Rb)}"))
        cmds.append(f"add line | x {pa[0]:.3f} y {pa[1]:.3f} z {pa[2]:.3f} "
                    f"| x {pb[0]:.3f} y {pb[1]:.3f} z {pb[2]:.3f} "
                    f"| color {color} | dashed false | type {lw} | radius {radius}")
    return cmds, used_sym, skipped, legend


# --------------------------------------------------------------------------- #
#  Main                                                                        #
# --------------------------------------------------------------------------- #

def main() -> None:
    ap = argparse.ArgumentParser(
        description="Layered base-pair notation + iCn3D commands to draw the pairs "
                    "in 3D, from a single annotated mmCIF.")
    ap.add_argument("cif", type=Path,
                    help="mmCIF with _ndb_base_pair_list/_annotation categories")
    ap.add_argument("chains", nargs="*",
                    help="chain ids in ruler order; default: all chains that pair")
    ap.add_argument("--name", default="RNA", help="label for the notation header")
    ap.add_argument("--id", default=None,
                    help="PDB id iCn3D should load (default: the cif file name); "
                         "must be an id iCn3D can load")
    ap.add_argument("--noncanonical", action="store_true",
                    help="draw only non-canonical pairs (cWW U-U/G-U included)")
    ap.add_argument("--block", type=int, nargs="?", const=BLOCK, default=None,
                    help=f"wrap notation lines into blocks of N columns "
                         f"(default {BLOCK})")
    ap.add_argument("--compact", action="store_true",
                    help="print sparse non-cWW layers as explicit pair lists "
                         "(e.g. A24,A31) instead of dot-bracket; saves space on large RNA")
    ap.add_argument("--script", metavar="FILE",
                    help="write the iCn3D commands to FILE (one per line) to load via "
                         "File > Open File > State/Script File, instead of a URL; "
                         "use this when there are too many pairs for a URL (~2000 char limit)")
    a = ap.parse_args()

    chains = a.chains or list_chains(a.cif)

    # 1) notation
    text, ok = build_notation(a.cif, chains, name=a.name, block=a.block,
                              compact=a.compact, noncanonical=a.noncanonical)
    print(text)
    print(f"\n# round-trip recovers all pairs exactly: {ok}", file=sys.stderr)

    # 2) 3D view
    cmds, used_sym, skipped, legend = build_icn3d_lines(a.cif, chains,
                                                        noncanonical=a.noncanonical)
    pdb_id = a.id or a.cif.stem.lower()   # pdbids load lowercase in iCn3D
    # 'set assembly on' must come first so the symmetry mate exists before the lines.
    draw_cmds = (["set assembly on"] if used_sym else []) + cmds

    if not cmds:
        msg = "non-canonical " if a.noncanonical else ""
        hint = " (drop --noncanonical to see the canonical pairs)" if a.noncanonical else ""
        print(f"\n# Nothing to draw: this structure has no {msg}pairs{hint}.",
              file=sys.stderr)
        sys.exit(0 if ok else 1)

    print(f"\n# iCn3D 3D view: {len(cmds)} pairs"
          f"{', non-canonical only' if a.noncanonical else ''}"
          f"{'; symmetry assembly included' if used_sym else ''}.")
    if legend:
        print("# non-canonical pairs (each drawn in its own colour):")
        for color, label in legend:
            print(f"#   {label}  #{color}")

    if a.script:
        # self-contained script file: load the structure, then draw the pairs.
        Path(a.script).write_text("\n".join([f"load pdb {pdb_id}"] + draw_cmds) + "\n")
        print(f"# Script written to {a.script} -> in iCn3D: File > Open File > State/Script File")
    else:
        print("#")
        print(f"# STEP 1 - open the structure in iCn3D (any browser):")
        print(f"#   {ICN3D_URL}?pdbid={pdb_id}")
        print("# STEP 2 - paste the lines below into iCn3D's command box (the '>' log at the")
        print("#   bottom), click right after the '>', and press Enter:")
        print("\n".join(draw_cmds))

    if skipped:
        print(f"# {skipped} pairs skipped (no coordinates)", file=sys.stderr)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
