# Layered base-pair notation

A lossless 2D dot-bracket notation for Watson-Crick **and** non-Watson-Crick base
pairs. It handles bases that pair more than once, pairs across chains, and
crystal-symmetry copies. Every example below round-trips: the notation parses
back to the exact FR3D pairs and Leontis-Westhof types.

## Method

- **Layers by LW type.** Each of the 18 directed LW types gets a fixed layer
  (L1 = cWW, L2 = cWH, ... L18 = tSS). Only the layers that actually occur are
  printed, and the layer label names the type, so there is no separate type list.
  A base that pairs twice with the same type spills to an overflow layer (L19+).
- **One ruler.** All residues are columns on a single axis, in order. A pair is
  one opening and one closing bracket on the same layer line. Crossing pairs
  within a layer use nested bracket levels `( [ { <`, the usual pseudoknot
  convention.
- **Multiple chains.** Chains are laid end to end on the ruler, separated by `&`
  (the ViennaRNA cofold strand separator). An inter-chain pair is just a bracket
  that opens in one chain's span and closes in another's.
- **Unambiguous residue identity.** Every residue is keyed by
  `(chain, symmetry, number)`. The symmetry field distinguishes a residue from a
  crystallographic copy of the same chain; the asymmetric-unit copy carries the
  explicit identity operator `1_555`. This makes self-complementary duplexes and
  crystal-packing contacts come out as separate strands rather than a residue
  pairing with itself.
- **Input sources.** Sequence and numbering come from `pdbx_poly_seq_scheme`,
  with a fallback to `_atom_site` (the coordinates) when that category is absent,
  e.g. modeling output. RNA and DNA are both handled; modified residues are shown
  in lowercase.
- **Scaling.** The number of layers is bounded by how many LW types occur (at
  most 18, plus a few overflow lines) and does not grow with sequence length;
  only the line length grows, like a FASTA line. A fixed-width mode (`--block N`)
  wraps long lines into blocks, FASTA/Stockholm style.

## Capabilities, phase by phase

### 1. Non-canonical pairs and multi-pairing (single chain)

9CFN, chain A. Several non-WC types are present, and 2 bases (5, 29) pair with
more than one partner. Each partner is kept on its own type layer instead of
being dropped, which is what a single dot-bracket line would have to do.

```
>9CFN complex: chains A[1_555](1-59)   ('&' separates chains)
seq         : AAGUACCCUCCAAGCCCUACAGGUUGGAAGAGGGGGCUAUCAGUCCUGUAGGCAGACUC
L1 cWW      : .(((.([[[[[..{.{{{{.{{{...)..].]]]])))......}}}.}}}}}......
L10 tWW     : .......................(......)............................
L11 tWH     : ........................(...)..............................
L12 tWS     : ....(....................)..(.............)................
L13 tHW     : ....(....................................).................
```

### 2. Pairs across chains

1XPE, the HIV-1 DIS kissing-loop dimer. Each chain folds into its own stem (the
`(((...)))`), and the two loops pair with each other through inter-chain pairs
(the `[[[ ... ]]]` crossing the `&`). They use the `[ ]` level because they cross
the intra-chain stems.

```
>1XPE complex: chains A[1_555](1-23), B[1_555](1-23)   ('&' separates chains)
seq         : CUUGCUGAAGCGCGCACGGCAAG&CUUGCUGAAGCGCGCACGGCAAG
L1 cWW      : (((((((..[[[[[[.)))))))&(((((((..]]]]]].)))))))
```

### 3. DNA

6NJQ, Structure of TBP-Hoogsteen containing DNA complex; chains C and D are a DNA duplex. 
Both Watson-Crick and a non-WC cWH pair occur between the chains, on separate type layers.

```
>6NJQ complex: chains C[1_555](201-214), D[1_555](215-228)   ('&' separates chains)
seq         : GCTATAAACGGGCA&TGCCCGTTTATAGC
L1 cWW      : (((.((((.((((.&.)))).)))).)))
L2 cWH      : ........(.....&.....)........
```

### 4. Crystal symmetry / self-complementary duplexes

2Q1R. The duplex partner is a symmetry copy: FR3D reports the pairs as chain A to
chain A, with a symmetry operator on one partner. Keying by
`(chain, symmetry, number)` renders this as two strands, `A[1_555]` and
`A[2_755]`, instead of a residue pairing with itself.

```
>2Q1R complex: chains A[1_555](1001-1012), A[2_755](1001-1012)   ('&' separates chains)
seq         : CGCGAAUUAGCG&CGCGAAUUAGCG
L1 cWW      : ((((((((((((&))))))))))))
```

### 5. Works without pdbx_poly_seq_scheme

The same 9CFN, built after removing the `pdbx_poly_seq_scheme` category, so the
sequence and numbering come from `_atom_site` alone. It round-trips to the same
pairs (residues 56-59 are unmodeled, have no coordinates, and are therefore
absent). This covers modeling results that lack that category.

```
>9CFN complex: chains A[1_555](1-55)   ('&' separates chains)
seq         : AAGUACCCUCCAAGCCCUACAGGUUGGAAGAGGGGGCUAUCAGUCCUGUAGGCAG
L1 cWW      : .(((.([[[[[..{.{{{{.{{{...)..].]]]])))......}}}.}}}}}..
L10 tWW     : .......................(......)........................
L11 tWH     : ........................(...)..........................
L12 tWS     : ....(....................)..(.............)............
L13 tHW     : ....(....................................).............
```

## Requirements

Python **3.10 or newer** (the script uses the `X | None` type-hint syntax, which
3.10 introduced); tested on Python 3.12. Only the Python standard library is
used, so there is nothing else to install. The easy way to get a recent Python
is a conda environment:

```
conda create -n lbn python=3.12
conda activate lbn
```

## Reproducing

```
python3 layered_basepairs.py <cif> <tsv> [chains...] [--name NAME] [--block N]
```

The chain ids are optional — with none given, all chains present in the FR3D TSV
are used automatically.

Base pairs come from the FR3D basepairs TSV (provided in `examples/`); the
matching mmCIF is downloaded from RCSB. For example, the 1XPE kissing-loop dimer:

```
wget https://files.rcsb.org/download/1XPE.cif
python3 layered_basepairs.py 1XPE.cif examples/1xpe_fr3d_basepairs.tsv A B --name 1XPE
# or let it pick up the chains itself:
python3 layered_basepairs.py 1XPE.cif examples/1xpe_fr3d_basepairs.tsv --name 1XPE
```

The notation goes to stdout; the round-trip check (`True` = lossless) is printed
to stderr.
