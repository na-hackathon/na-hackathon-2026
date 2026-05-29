"""Layered base-pair notation from an FR3D base-pair TSV.

Chains are laid end to end on one ruler, separated by '&' (the ViennaRNA cofold
convention), and a pair between two chains is a bracket that opens in one chain
and closes in another. Each base is keyed by (chain, symmetry, number): the
symmetry field distinguishes a residue from a symmetry-generated copy of the
same chain, so self-complementary duplexes and crystal-packing contacts (where
FR3D reports a pair as chain A to chain A with a symmetry operator on one
partner) are represented as two strands rather than a spurious self-pairing.

The sequence and residue numbers come from pdbx_poly_seq_scheme when present,
otherwise from _atom_site (the coordinates), so structures lacking that
category (e.g. modelling results) still work.

This file owns the FR3D-specific reading (unit-id parsing, pair list); the rest
-- constants, sequence/numbering reading, the layered-notation builder, parser
helpers -- lives in common.py.

Usage:
    python3 layered_basepairs.py <cif> <tsv> [chains...] [options]

The mmCIF is for sequence + residue numbers only (any vanilla CIF works).
The FR3D TSV provides the base-pair list and Leontis-Westhof families.

Examples -- most users only need the first one:

    # 1. Simplest. --name auto-detected from the CIF; chains auto-picked from
    #    every chain that appears in the TSV.
    python3 layered_basepairs.py 9CFN.cif 9cfn_fr3d_basepairs.tsv

    # 2. Also list residues that have no pair (always derived from the TSV
    #    here -- FR3D has no separate 'unpaired list' source).
    python3 layered_basepairs.py 9CFN.cif 9cfn_fr3d_basepairs.tsv --unpaired

    # 3. Two-chain complex, pick the chain order on the ruler.
    python3 layered_basepairs.py 1XPE.cif 1xpe_fr3d_basepairs.tsv A B

    # 4. Non-canonical pairs only, with slot-numbered row labels (L0, L1, ...).
    python3 layered_basepairs.py 9CFN.cif 9cfn_fr3d_basepairs.tsv --noncanonical --layer

    # 5. Long sequences: wrap into 100-column blocks AND keep non-canonical
    #    rows as a short pair list instead of full-width dot-bracket.
    python3 layered_basepairs.py 9CFN.cif 9cfn_fr3d_basepairs.tsv --block --compact

    # 6. Add a '# chains: ...' header comment describing each strand.
    python3 layered_basepairs.py 2Q1R.cif 2q1r_fr3d_basepairs.tsv --metadata

    # 7. Override the auto-detected name.
    python3 layered_basepairs.py 9CFN.cif 9cfn_fr3d_basepairs.tsv --name MY-LABEL

The notation prints to stdout; the round-trip result (True = lossless) to stderr.
"""

import argparse
import sys
from pathlib import Path

from common import (
    BLOCK, IDENTITY,
    _flip_lw, build_notation, read_entry_id, read_residues,
)


# --------------------------------------------------------------------------- #
#  FR3D-specific reading                                                       #
# --------------------------------------------------------------------------- #

def _sym(unit: list[str]) -> str:
    """Symmetry operator from an FR3D unit id. A missing operator denotes the
    asymmetric-unit copy, which is the crystallographic identity, so it is
    reported explicitly as 1_555. FR3D unit id format is
    pdb|model|chain|comp|num|atom|alt|ins|symmetry."""
    return unit[8] if len(unit) > 8 and unit[8] else IDENTITY


def list_chains(tsv: Path) -> list[str]:
    """Chains that appear in the FR3D pair list, sorted. Used as the default
    when no chains are given on the command line."""
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


def read_pairs(tsv: Path, chains: set[str]) -> list[tuple[tuple, tuple, str]]:
    """De-duplicated pairs (Ra, Rb, lw), where each residue R is
    (chain, symmetry, number). Both ends must be among `chains`; ordering is
    canonical so the LW direction is consistent."""
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


# --------------------------------------------------------------------------- #
#  Main                                                                        #
# --------------------------------------------------------------------------- #

if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Layered base-pair notation for an RNA complex "
                    "(symmetry-aware), driven by an FR3D base-pair TSV.")

    # Positional --------------------------------------------------------------
    ap.add_argument("cif", type=Path,
                    help="input mmCIF (used only for sequence + residue numbers)")
    ap.add_argument("tsv", type=Path, help="FR3D base-pair TSV file")
    ap.add_argument("chains", nargs="*",
                    help="chain IDs in ruler order; "
                         "default: every chain that appears in the TSV")

    # Identity ----------------------------------------------------------------
    ap.add_argument("--name", default=None,
                    help="header label shown in the '>...' line; "
                         "default: read from the CIF itself (_entry.id, or "
                         "the 'data_XXXX' block name) -- so for 8BWT.cif you "
                         "do not have to type --name 8BWT")

    # Display filters ---------------------------------------------------------
    ap.add_argument("--noncanonical", action="store_true",
                    help="drop canonical Watson-Crick pairs (A-U, G-C, A-T) "
                         "from the display; non-canonical cWW pairs "
                         "(U-U / U-G wobbles) stay on the cWW row")
    ap.add_argument("--compact", action="store_true",
                    help="keep the canonical WC layer as full-width dot-bracket, "
                         "but print every non-canonical layer as a short pair "
                         "list (e.g. 'A24,A31  A25,A29'); saves vertical space "
                         "on large RNAs")
    ap.add_argument("--block", type=int, nargs="?", const=BLOCK, default=None,
                    help=f"wrap each row into blocks of N columns "
                         f"(default {BLOCK} when --block is given with no value); "
                         f"useful for very long sequences")

    # Row / strand labelling --------------------------------------------------
    ap.add_argument("--layer", action="store_true",
                    help="prepend slot numbers to each row label "
                         "(L0 WC, L1 cWW, L2 cWH, ... L18 tSS); "
                         "off by default -- labels are plain family names")

    # Extra comment lines below the header -----------------------------------
    ap.add_argument("--metadata", action="store_true",
                    help="add a '# chains: ...' line describing each strand "
                         "(chain, symmetry operator, residue range)")
    ap.add_argument("--unpaired", action="store_true",
                    help="add a '# unpaired (N): ...' line listing residues "
                         "with no pair. Always derived from the FR3D TSV "
                         "here -- FR3D has no separate 'unpaired list' source")
    a = ap.parse_args()

    chains = a.chains or list_chains(a.tsv)
    per_chain = read_residues(a.cif, chains)
    pairs = read_pairs(a.tsv, set(chains))
    # If --name wasn't given, pull it from the CIF so users don't have to retype
    # the PDB id that the file already declares about itself.
    name = a.name or read_entry_id(a.cif) or "RNA"
    text, ok = build_notation(per_chain, pairs, chains, name=name,
                              block=a.block, compact=a.compact,
                              noncanonical=a.noncanonical,
                              show_layer=a.layer, show_metadata=a.metadata,
                              show_unpaired=a.unpaired,
                              unpaired_source=("derived from FR3D pair list"
                                               if a.unpaired else ""))
    print(text)
    print(f"\n# round-trip recovers all pairs exactly: {ok}", file=sys.stderr)
    sys.exit(0 if ok else 1)
