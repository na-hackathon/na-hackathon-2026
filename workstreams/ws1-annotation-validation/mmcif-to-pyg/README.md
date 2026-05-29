# mmcif-to-pyg — RNA base-pair graphs & GNN models

Turn RNA base-pair data (mmCIF) into graphs and learn on them. Two input formats are
supported, both producing the same **rnaglib-style 2.5D graph** (nodes = nucleotides,
edges = backbone + Leontis–Westhof base pairs), which then feeds Relational-GCN models
for three tasks: base-pair **typing**, **contact** prediction, and joint **detect + type**.

Everything runs in a reproducible [pixi](https://pixi.sh) environment (Python 3.12; torch,
torch_geometric, gemmi, rnaglib, fr3d, forgi, scipy, seaborn, jupyterlab).

---

## Setup

```bash
cd workstreams/ws1-annotation-validation/mmcif-to-pyg
pixi install                      # creates .pixi/envs/default
pixi run python -c "import torch, torch_geometric, rnaglib, gemmi; print('ok')"
```

Run anything with `pixi run <cmd>` (it activates the env). For notebooks, either
`pixi run jupyter lab` or select the registered kernel **"Python (mmcif-to-pyg pixi)"**.

---

## Data formats

| format | example | has coordinates? | how it's read |
|---|---|---|---|
| **NDB / NAPAIR** annotation mmCIF | `data/00008bwt.mmcif` | no (`_ndb_base_pair_*` only) | `mmcif_to_pyg.py`, `ndb_to_rnaglib.py` |
| **DNATCO "extended"** mmCIF | fetched from dnatco.datmos.org | yes (`_atom_site` + NDB annotations) | `dnatco_to_graph.py` |
| **rnaglib graph JSON** | produced by the adapters | per-nucleotide `xyz_C1p/glyN/P` | `rnaglib.utils.load_graph` |

A graph node carries `nt_code` (A/C/G/U/N), `chain_id`, `index`, `is_modified`, and (for
DNATCO) `xyz_C1p` / `xyz_glyN` / `xyz_P`. Edges carry an `LW` label: backbone `B53`/`B35`
or a base-pair family (`cWW`, `tWS`, …).

---

## Pipeline (DNATCO → graph → train)

```bash
# 1. fetch DNATCO coordinate mmCIFs  (ids from RCSB, a file, or the CLI)
pixi run python dnatco_fetch.py --rcsb-rna 1000 -o cifs/
#    or from a list: cut -d, -f1 ../../../data/RS25_RNA_c2.csv | sort -u > ids.txt
#    pixi run python dnatco_fetch.py --ids-file ids.txt -o cifs/

# 2. convert to rnaglib graph JSON (coordinates included)
pixi run python dnatco_to_graph.py cifs/ -o graphs/

# 3. train a model (examples)
pixi run python train_pair_detect_type.py graphs/ --epochs 25      # joint detect + type
pixi run python train_rs25.py graphs/ --epochs 30                  # detect+type, train/val/test
```

The annotation-only path (no fetch needed) for the committed example:

```bash
pixi run python mmcif_to_pyg.py   ../../../data/00008bwt.mmcif -o graphs_8bwt.pt   # -> PyG graphs
pixi run python ndb_to_rnaglib.py ../../../data/00008bwt.mmcif -o 8bwt.json        # -> rnaglib graph
```

---

## Files

### Library
- **`rna_pairs.py`** — shared core: constants, `PairData`, `parse_graph`, `geom_features`,
  `sample_hard_negatives`, the configurable `PairGNN`, training loops, metrics, and
  `is_canonical` (base-aware canonical test). The trainers are thin wrappers over this.

### Data → graph
- **`mmcif_to_pyg.py`** — NDB annotation mmCIF → PyTorch Geometric `Data` (one graph per PDB model).
- **`ndb_to_rnaglib.py`** — NDB annotation mmCIF → rnaglib-loadable 2.5D graph (so rnaglib's
  own `load_graph` / drawing work on coordinate-free annotation files).
- **`dnatco_fetch.py`** — download DNATCO extended mmCIFs (concurrent; `--rcsb-rna N` /
  `--ids-file` / CLI ids).
- **`dnatco_to_graph.py`** — DNATCO mmCIF → rnaglib graph JSON with backbone + LW base-pair
  edges and representative coordinates (C1′, glycosidic N, P).

### Models / training
| script | task | input | output |
|---|---|---|---|
| `train_pair_interactions.py` | **typing** of given pairs | identity + backbone + a generic `PAIR` edge | LW family (18-way) |
| `train_contact_prediction.py` | **contact** (which residues pair) | identity + backbone + geometry; hard negatives | pair / no-pair (`--no-geom` ablation) |
| `train_pair_detect_type.py` | **detect + type** (joint) | identity + backbone + geometry | NONE + 18 LW families |
| `train_rs25.py` | detect+type with **train/val/test** | a graph dir | best-val checkpoint, test metrics |
| `train_contact_seqonly.py` | **coordinate-free** contact | sequence + backbone only | pair / no-pair (`--with-geom` ablation) |
| `train_dnatco.py` | self-supervised masked-nucleotide | identity (masked) + LW + backbone | recover A/C/G/U |

### Notebooks
- **`dnatco_pipeline_demo.ipynb`** — end-to-end on the RS25 set: fetch → graph → the three
  tasks, with a structure-level train/val/test split and per-class / per-structure plots.
- **`walkthrough/pipeline_walkthrough.ipynb`** — visualizes each pipeline stage (mmCIF → graph
  → arc diagram → backbone layout → 3D coordinates → candidate-pair geometry).
- **`mmcif_to_pyg.ipynb`**, **`ndb_via_rnaglib.ipynb`**, **`rnaglib_8bwt.ipynb`** — focused
  demos of the conversion / rnaglib paths.

### Tests
```bash
pixi run python -m unittest test_mmcif_to_json -v
```
Covers mmCIF→JSON conversion + round-trip, which information is preserved vs dropped, and
building/visualizing the graph with networkx.

---

## Leontis–Westhof families & "canonical"

Each base pair has a 3-char LW code = **orientation** (`c`is/`t`rans) + the two **edges**
used (`W` Watson-Crick, `H` Hoogsteen, `S` Sugar) → 12 geometric families (18 directed
labels in the data).

**Canonical** here is base-aware (`rna_pairs.is_canonical`): a pair is canonical only if it
is **cWW *and* a Watson-Crick base combination** (A·U / G·C; G·U wobble optional). A `cWW_U-U`
is therefore **non-canonical** — cWW geometry but not a Watson-Crick pair. Everything else
(other families, or cWW with non-WC bases) is non-canonical — the WS2-relevant signal.

---

## Notes & caveats

- **CPU only** here (no GPU); training thousands of large graphs (e.g. ribosomes) is slow.
- **Geometry tasks re-learn FR3D.** Contact / detect+type are *given* the 3D geometry, while
  the labels are derived from that geometry — so high typing accuracy is expected; it's
  structure *annotation*, not de-novo prediction. `train_contact_seqonly.py` is the honest
  coordinate-free (prediction) setting.
- **Detection metrics use a curated candidate set** (positives + sampled hard negatives), not
  all O(N²) pairs — read F1 accordingly.
- **Random structure splits can leak homologs** (RS25 has many related RNAs); a
  sequence/family-clustered split gives a less optimistic estimate.

Generated artifacts (`.pixi/`, fetched `cifs/`, `graphs/`, `demo/`, `rs25/`, `*.pt`,
`*.rnaglib.json`, logs) are gitignored.
