# Layered base-pair notation

A lossless 2D dot-bracket notation for Watson-Crick **and** non-Watson-Crick base
pairs. It handles bases that pair more than once, pairs across chains, and
crystal-symmetry copies. Every example below round-trips: the notation parses
back to the exact pairs and Leontis-Westhof types.

## Method

- **Layers by LW type, canonical on its own row.** The canonical Watson-Crick
  pairs (cWW geometry AND A-U / G-C / A-T bases) sit on a dedicated `WC` row.
  Non-canonical cWW pairs (cWW geometry but non-WC bases, e.g. U-U or U-G
  wobbles) sit on a separate `cWW` row. The remaining 17 directed LW families
  each get their own row (`tWW`, `cWH`, `tWH`, ..., `tSS`). Only the rows that
  actually occur are printed, and the row label names the family. A base that
  pairs twice with the same family spills to an overflow row. With `--layer`
  the rows additionally carry slot numbers (`L0 WC`, `L1 cWW`, `L10 tWW`, ...).
- **One ruler.** All residues are columns on a single axis, in order. A pair is
  one opening and one closing bracket on the same row. Crossing pairs within a
  row use nested bracket levels `( [ { <`, the usual pseudoknot convention.
- **Multiple chains.** Chains are laid end to end on the ruler, separated by
  `&` (the ViennaRNA cofold strand separator). An inter-chain pair is a bracket
  that opens in one chain's span and closes in another's.
- **Unambiguous residue identity.** Every residue is keyed by
  `(chain, symmetry, number)`. The symmetry field distinguishes a residue from
  a crystallographic copy of the same chain; the asymmetric-unit copy carries
  the explicit identity operator `1_555`. This makes self-complementary
  duplexes and crystal-packing contacts come out as separate strands rather
  than a residue pairing with itself.
- **Parseable header.** NCBI/FASTA-style: `>name|<strand>|<strand>|...`. Each
  strand is a chain id; symmetry mates append `:operator`
  (e.g. `>2Q1R|A|A:2_755`). Identity `1_555` is implicit. Use `--metadata` to
  add a `# chains: A[1_555](1-59); ...` comment line below the header.
- **Input sources.** Sequence and numbering come from `pdbx_poly_seq_scheme`,
  with a fallback to `_atom_site` (the coordinates) when that category is
  absent, e.g. modelling output. RNA and DNA are both handled; modified
  residues are shown in lowercase.
- **Scaling.** The number of rows is bounded by how many LW families occur (at
  most 1 WC + 18 LW + a few overflow lines) and does not grow with sequence
  length; only the line length grows, like a FASTA line. A fixed-width mode
  (`--block N`) wraps long lines into blocks, FASTA/Stockholm style.

## Files

- `common.py` — shared constants, CIF-sequence reading, layout / render / parse
  helpers, and the top-level `build_notation` builder. Both scripts below
  import from it.
- `layered_basepairs.py` — notation from an **FR3D base-pair TSV** plus an
  mmCIF (used for sequence and residue numbers). Symmetry-aware via FR3D
  unit-id parsing.
- `standalone_lbn_script.py` — notation directly from a **DNATCO-extended
  mmCIF** (`_ndb_base_pair_list` + `_ndb_base_pair_annotation`), plus an iCn3D
  3D overlay so you can see every pair on the real structure. No external pair
  list needed.
- `notation_to_cif.py` — the **inverse**: takes a layered-notation file plus
  any mmCIF of the same structure and writes a new mmCIF whose base-pair
  annotation is rebuilt straight from the notation, proving the round-trip
  in both directions.

## Capabilities, phase by phase

### 1. Non-canonical pairs and multi-pairing (single chain)

9CFN, chain A. Several non-WC types are present, and 2 bases (5, 29) pair with
more than one partner. Each partner is kept on its own family row instead of
being dropped, which is what a single dot-bracket line would have to do. Note
the U-G wobbles at A4-A36 and A9-A33 sit on `cWW` (non-canonical), not on `WC`.

```
>9CFN|A
seq         : AAGUACCCUCCAAGCCCUACAGGUUGGAAGAGGGGGCUAUCAGUCCUGUAGGCAGACUC
WC          : .((..([[.[[..{.{{{{.{{{...)..].].]].))......}}}.}}}}}......
cWW         : ...(....(.......................)..).......................
tWW         : .......................(......)............................
tWH         : ........................(...)..............................
tWS         : ....(....................)..(.............)................
tHW         : ....(....................................).................
```

### 2. Pairs across chains

1XPE, the HIV-1 DIS kissing-loop dimer. Each chain folds into its own stem
(the `(((...)))`), and the two loops pair with each other through inter-chain
pairs (the `[[[ ... ]]]` crossing the `&`). They use the `[ ]` level because
they cross the intra-chain stems. The non-canonical row picks up one U-G
wobble per chain.

```
>1XPE|A|B
seq         : CUUGCUGAAGCGCGCACGGCAAG&CUUGCUGAAGCGCGCACGGCAAG
WC          : (((((.(..[[[[[[.).)))))&(((((.(..]]]]]].).)))))
cWW         : .....(...........).....&.....(...........).....
```

### 3. DNA

6NJQ, structure of a TBP-Hoogsteen-containing DNA complex; chains C and D are
a DNA duplex. Both Watson-Crick and a non-WC `cWH` pair occur between the
chains, on separate family rows.

```
>6NJQ|C|D
seq         : GCTATAAACGGGCA&TGCCCGTTTATAGC
WC          : (((.((((.((((.&.)))).)))).)))
cWH         : ........(.....&.....)........
```

### 4. Crystal symmetry / self-complementary duplexes

2Q1R. The duplex partner is a symmetry copy: FR3D reports the pairs as chain A
to chain A, with a symmetry operator on one partner. Keying by
`(chain, symmetry, number)` renders this as two strands, the identity copy `A`
and the symmetry copy `A:2_755`, instead of a residue pairing with itself.

```
>2Q1R|A|A:2_755
seq         : CGCGAAUUAGCG&CGCGAAUUAGCG
WC          : (((.((((.(((&))).)))).)))
cWW         : ...(....(...&...)....)...
```

### 5. Works without `pdbx_poly_seq_scheme`

The same 9CFN, built after removing the `pdbx_poly_seq_scheme` category, so
the sequence and numbering come from `_atom_site` alone. It round-trips to the
same pairs (residues 56-59 are unmodeled, have no coordinates, and are
therefore absent). This covers modelling results that lack that category.

```
>9CFN|A
seq         : AAGUACCCUCCAAGCCCUACAGGUUGGAAGAGGGGGCUAUCAGUCCUGUAGGCAG
WC          : .((..([[.[[..{.{{{{.{{{...)..].].]].))......}}}.}}}}}..
cWW         : ...(....(.......................)..)..................
tWW         : .......................(......)........................
tWH         : ........................(...)..........................
tWS         : ....(....................)..(.............)............
tHW         : ....(....................................).............
```

### 6. Unpaired residues

A `.` in a bracket row means "this residue has no partner **on that layer**"
— that is a per-row statement. With `--unpaired`, the script adds a single
`#`-prefixed comment line summarising "this residue has no partner **anywhere
in the structure**". Both views coexist; one is per-layer, the other is
structural.

```
$ python3 standalone_lbn_script.py 9HRF.cif --unpaired
>9HRF|A
# unpaired (13): A8, A11, A29-32, A44, A49-51, A61, A69-70  [derived from pair list]
seq         : GGGCGCAUUUUGGUAGGUCGGUCGCUGCUUCGGCAGUGAGGGGUAGGCAUUGCUGGCCUAGGGUGCCCUU
WC          : ((((.(...(...[[[[.[.[((.((((....)))).)).....)..(...)].].]]]]..).))))..
...
```

Consecutive residue numbers within one strand are collapsed to ranges
(`A29-32` instead of `A29, A30, A31, A32`); symmetry mates use the same
bracketed-operator form the rest of the script uses (e.g. `A2[2_565]`).

**Source rule.** The script picks where to get the unpaired list using a
simple two-step check:

1. **If the CIF contains `_ndb_base_unpaired_list`** — DNATCO's own
   per-model annotation — read it directly (model 1 by default). This is the
   authoritative source, particularly for NMR ensembles where the unpaired
   set varies across models.
2. **Otherwise** — derive it from `_ndb_base_pair_list`: every residue in
   the chain that never appears in any pair row is unpaired.

The output records which source was used:

```
# unpaired (4): A12-15  [from _ndb_base_unpaired_list, model 1]
# unpaired (2): A13-14  [derived from pair list]
```

The bundled `examples/8VJT.cif` is the only example file shipped with the
`_ndb_base_unpaired_list` block already populated, so it is the easiest way
to try the "from `_ndb_base_unpaired_list`" path:

```
python3 standalone_lbn_script.py examples/8VJT.cif --unpaired
# >8VJT|A|A:2_565|...
# # unpaired (2): A2, B2  [from _ndb_base_unpaired_list, model 1]
```

The same 8VJT through the FR3D path (which has no separate unpaired list)
derives the unpaired set from the TSV instead, and reports it across every
displayed strand including the symmetry mates:

```
python3 layered_basepairs.py examples/8VJT.cif examples/8vjt_fr3d_basepairs.tsv --unpaired
# # unpaired (8): A2, A2[2_565], A2[3_555], A2[4_455], B2, ...  [derived from FR3D pair list]
```

Both views are correct for the scope they describe (asymmetric unit vs the
displayed assembly).

The FR3D-driven script (`layered_basepairs.py`) always derives, since FR3D
ships only a pair list.

**Why a `#` comment line and not a new bracket row?** Standard 2D viewers
(VARNA, R2DT, forna, fornac, ViennaRNA) require `.` for unpaired in
dot-bracket. Putting a custom marker on a bracket row would either break
bracket balance (`(` for unpaired) or break tool compatibility (`_`, `-`).
A `#` comment follows the FASTA / Stockholm convention, is ignored by every
standard tool, and stays human-readable. It also round-trips correctly:
`notation_to_cif.py` ignores `#` lines, so the comment never affects the
rebuilt CIF.

The display-filter is honoured **but the count is structural**, by design:
running `--unpaired --noncanonical` reports the same residues as `--unpaired`
alone, because a residue paired only via Watson-Crick is still paired in the
structure — its pair is just hidden when the WC layer is dropped from view.

## `standalone_lbn_script.py` — also draws the pairs in 3D

If your CIF already carries the DNATCO base-pair annotation
(`_ndb_base_pair_list` + `_ndb_base_pair_annotation`), this script reads the
pairs straight from it (no external TSV needed) and, in the same run, prints
iCn3D `add line` commands that draw every pair on the real 3D structure.
Canonical WC pairs are drawn thin grey; each non-canonical pair gets its own
bold colour (listed in the legend below the notation).

```
python3 standalone_lbn_script.py 9cfn_dnatco.cif A --name 9CFN
```

Sample output (truncated):

```
>9CFN|A
seq         : AAGUACCCUCCAAGCCCUACAGGUUGGAAGAGGGGGCUAUCAGUCCUGUAGGCAGACUC
WC          : .((..([[.[[..{.{{{{.{{{...)..].].]].))......}}}.}}}}}......
cWW         : ...(....(.......................)..).......................
tWW         : .......................(......)............................
tWH         : ........................(...)..............................
tWS         : ....(....................)..(.............)................
tHW         : ....(....................................).................

# iCn3D 3D view: 22 pairs.
# non-canonical pairs (each drawn in its own colour):
#   cWW U-G A4-A36  #FF00FF
#   tWS A-G A5-A26  #FF8000
#   tHW A-A A5-A42  #00C800
#   cWW U-G A9-A33  #00E5E5
#   tWW U-A A24-A31  #FFD000
#   tWH U-A A25-A29  #FF0000
#   tWS A-G A29-A43  #0000FF
#
# STEP 1 - open the structure in iCn3D (any browser):
#   https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?pdbid=9cfn
# STEP 2 - paste the lines below into iCn3D's command box (the '>' log at the
#   bottom), click right after the '>', and press Enter:
add line | x ... | x ... | color FF00FF | dashed false | type cWW | radius 0.8
add line | x ... | x ... | color 888888 | dashed false | type cWW | radius 0.3
...
```

For structures with too many pairs to paste comfortably, use `--script
9cfn.icn3d` to write a script file and load it via iCn3D's
**File > Open File > State/Script File**. Symmetry-mate structures get an
automatic `set assembly on` directive so the operator copy exists before the
lines are drawn.

## `notation_to_cif.py` — the inverse direction

The forward scripts prove notation is recoverable from a CIF. This script
proves a CIF is recoverable from notation: together they close the loop and
make the layered notation a true lossless representation of base-pair
annotation. Whatever you can write down in dot-bracket layers can be put back
into an mmCIF and read by any tool that understands mmCIF, with the same pair
set.

```
python3 notation_to_cif.py <notation.txt> <any.cif> -o <new.cif>
```

The script regenerates the three categories that carry the base-pair
information in DNATCO-extended mmCIF:

- `_ndb_base_pair_list` — who pairs with whom (auth numbering + operator id)
- `_ndb_base_pair_annotation` — the Leontis-Westhof family per pair
- `_pdbx_struct_oper_list` — maps each operator id to its `1_555`-style name

The two annotation loops are rebuilt from the notation. The operator list is
left untouched if the input CIF already lists every operator the notation
needs (matrix data preserved), or extended when the notation references an
operator the CIF does not have (e.g. a `2_755` symmetry mate applied to a
vanilla CIF).

The input CIF does **not** have to be DNATCO-extended. The script falls back
through three metadata sources -- existing `_ndb_base_pair_list` rows, then
`_pdbx_poly_seq_scheme`, then `_atom_site` -- so a vanilla PDB-downloaded CIF
or a minimal atom-site-only CIF works too.

### Round-trip verification

After writing the new CIF the script re-reads it with the standard reader and
confirms the pair set `(chain, num, chain, num, LW)` matches what the notation
encoded:

```
$ python3 standalone_lbn_script.py 9cfn_dnatco.cif A --name 9CFN > 9CFN.lbn
$ python3 notation_to_cif.py 9CFN.lbn 9cfn_dnatco.cif -o 9CFN_rebuilt.cif
# recovered 22 pairs from notation
# wrote 9CFN_rebuilt.cif
# notation -> CIF -> pairs round-trip: True  (22 pairs)
# rebuilt CIF == original CIF (semantically; all pairs preserved).
```

Verified on the example structures (8BWT, 1XPE, 6NJQ, 2Q1R, 9CFN, and
9CFN-without-`pdbx_poly_seq_scheme`) in default / `--compact` / `--layer` /
`--metadata` / `--noncanonical` modes — all True. Even stronger: re-running
the forward script on the rebuilt CIF yields a notation that is
byte-identical to the original.

### Header and label compatibility

The parser accepts both the new NCBI-style header (`>9CFN|A`,
`>2Q1R|A|A:2_755`) and the older verbose form
(`>9CFN complex: chains A[1_555](1-59) ...`). It also accepts both the default
family-only labels (`WC`, `cWW`, `tWW`, ...) and the `--layer` slot-numbered
ones (`L0 WC`, `L1 cWW`, `L10 tWW`, ...). The `# chains: ...` metadata comment
line is ignored.

### Limitation

The notation only encodes the pair list `(chain, num, chain, num, LW)`. It
does not encode 3D coordinates or operator matrices, so the rebuilt CIF
preserves whatever atoms and matrices were already in the input CIF and
replaces only the base-pair annotation blocks. Rebuilding atoms from scratch
would need structural modelling and is out of scope here.

## Common flags

Both forward scripts share these flags:

| Flag | Effect |
| --- | --- |
| `--compact` | Keep only the canonical `WC` row as full dot-bracket; print every non-canonical row as an explicit pair list (e.g. `A24,A31`). Saves space on large RNA. |
| `--noncanonical` | Drop true Watson-Crick pairs entirely (the `WC` row becomes absent). cWW U-U / U-G wobbles still appear on the `cWW` row. |
| `--layer` | Prepend slot numbers to each row label (`L0 WC`, `L1 cWW`, `L10 tWW`, ...). Off by default. |
| `--metadata` | Add a `# chains: A[1_555](1-59); '&' separates chains` comment line below the header. The line starts with `#` so any parser ignores it. |
| `--unpaired` | Add a `# unpaired (N): A8, A11, A29-32, ...` comment line listing residues that have no pair. Consecutive residues within a strand are collapsed to ranges. See **Unpaired residues** below for how the source is chosen. |
| `--block N` | Wrap lines into fixed-width blocks of N columns (handy for long sequences). `--block` with no value uses 100 columns. |
| `--name NAME` | Header label. **Default: read from the CIF itself** (`_entry.id`, falling back to the `data_XXXX` block name at the top of the file). So for `8BWT.cif` you get `>8BWT|...` without typing `--name`. |

## Requirements

Python **3.10 or newer** (the script uses the `X | None` type-hint syntax,
which 3.10 introduced); tested on Python 3.12. Only the Python standard
library is used, so there is nothing else to install. The easy way to get a
recent Python is a conda environment:

```
conda create -n lbn python=3.12
conda activate lbn
```

## Reproducing

### From an FR3D TSV (any mmCIF for sequence/numbering)

```
python3 layered_basepairs.py <cif> <tsv> [chains...] [options]
```

Optional flags: `--name`, `--block`, `--compact`, `--noncanonical`, `--layer`,
`--metadata`, `--unpaired`. With no chains, every chain present in the FR3D
TSV is used automatically. With no `--name`, the label is read from the
CIF's `_entry.id`.

Base pairs come from the FR3D basepairs TSV (provided in `examples/`); the
matching mmCIF is downloaded from RCSB. For example, the 1XPE kissing-loop
dimer:

```
wget https://files.rcsb.org/download/1XPE.cif
python3 layered_basepairs.py 1XPE.cif examples/1xpe_fr3d_basepairs.tsv A B
# or let it pick up the chains itself:
python3 layered_basepairs.py 1XPE.cif examples/1xpe_fr3d_basepairs.tsv
```

### From a DNATCO-extended mmCIF (no external pair list)

```
python3 standalone_lbn_script.py <cif> [chains...] [options]
```

Optional flags: `--name`, `--id`, `--compact`, `--noncanonical`, `--layer`,
`--metadata`, `--unpaired`, `--block`, `--script FILE`. With no chains,
every chain that appears in the pair list is used.

- `--name` defaults to `_entry.id` from the CIF (or the `data_XXXX` block name).
- `--id` overrides the PDB id iCn3D should load (defaults to the CIF file stem,
  lowercased).
- `--unpaired` reads `_ndb_base_unpaired_list` when the CIF contains it
  (DNATCO's own per-model annotation), otherwise derives from the pair list.

Minimal invocation, everything auto-detected:

```
python3 standalone_lbn_script.py 8BWT.cif --unpaired
```

The notation goes to stdout; the round-trip check (`True` = lossless) is
printed to stderr. The standalone additionally emits iCn3D `add line` commands
(or, with `--script FILE`, writes them to a self-contained iCn3D script).
