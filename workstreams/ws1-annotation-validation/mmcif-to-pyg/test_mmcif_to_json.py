#!/usr/bin/env python3
"""
Unit tests for the mmCIF -> rnaglib JSON graph conversion.

Verifies three things the project cares about:
  1. we can convert an mmCIF to a JSON graph (and rnaglib can read it back);
  2. which information is *used* (preserved) vs *lost* (dropped) in the JSON;
  3. the graph builds into a networkx graph and can be visualized.

Run:
    pixi run python -m unittest test_mmcif_to_json -v
"""
from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import matplotlib
matplotlib.use("Agg")          # headless: no display needed for the viz test
import matplotlib.pyplot as plt
import networkx as nx

from rnaglib.utils import dump_json, load_graph
from ndb_to_rnaglib import ndb_mmcif_to_rnaglib_graph

REPO = Path(__file__).resolve().parents[3]
NDB_MMCIF = REPO / "data" / "00008bwt.mmcif"          # NDB annotation-only mmCIF (committed)
SEQ_8BWT = "GGCGUUUUCGCUUCGGCGUUUACGCC"               # 8BWT chain A, 5'->3'


def _find_dnatco_cif() -> Path | None:
    """A locally cached DNATCO coordinate cif, if any (for the optional coord test)."""
    here = Path(__file__).resolve().parent
    for d in (here / "rs25" / "cifs", here / "demo" / "cifs"):
        if d.is_dir():
            for f in sorted(d.glob("*.cif")):
                return f
    return None


class TestNdbMmcifToJson(unittest.TestCase):
    """mmCIF (annotation-only) -> rnaglib JSON graph."""

    @classmethod
    def setUpClass(cls):
        assert NDB_MMCIF.exists(), f"missing test input {NDB_MMCIF}"
        cls.tmp = tempfile.TemporaryDirectory()
        cls.G = ndb_mmcif_to_rnaglib_graph(NDB_MMCIF, model=1)
        cls.json_path = Path(cls.tmp.name) / "8bwt.json"
        dump_json(str(cls.json_path), cls.G)
        cls.reloaded = load_graph(str(cls.json_path))
        cls.raw = NDB_MMCIF.read_text()
        cls.js = cls.json_path.read_text()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    # --- 1. conversion + round-trip -----------------------------------------
    def test_json_created_and_loadable(self):
        self.assertTrue(self.json_path.exists())
        self.assertGreater(self.json_path.stat().st_size, 0)
        self.assertIsInstance(self.reloaded, nx.DiGraph)
        self.assertEqual(self.reloaded.graph["pdbid"], "8bwt")

    def test_roundtrip_preserves_counts(self):
        self.assertEqual(self.reloaded.number_of_nodes(), self.G.number_of_nodes())
        self.assertEqual(self.reloaded.number_of_edges(), self.G.number_of_edges())

    # --- 2a. information USED (preserved) -----------------------------------
    def test_sequence_and_node_attrs_preserved(self):
        G = self.reloaded
        self.assertEqual(G.number_of_nodes(), 26)
        for _, d in G.nodes(data=True):
            for key in ("nt_code", "chain_id", "index"):
                self.assertIn(key, d)
        seq = "".join(G.nodes[n]["nt_code"] for n in sorted(G.nodes, key=lambda x: G.nodes[x]["index"]))
        self.assertEqual(seq, SEQ_8BWT)

    def test_base_pairs_and_backbone_preserved(self):
        G = self.reloaded
        lw = [d["LW"] for *_, d in G.edges(data=True)]
        self.assertIn("cWW", lw)                          # LW family kept
        self.assertIn("B53", lw)
        self.assertIn("B35", lw)                          # backbone kept
        n_pairs = sum(l not in ("B53", "B35") for l in lw) // 2
        self.assertEqual(n_pairs, 11)                     # model 1 has 11 base pairs

    # --- 2b. information LOST (dropped by the conversion) --------------------
    def test_validation_and_edge_detail_dropped(self):
        # the source mmCIF DOES contain these fields ...
        for present in ("napair_rmsd", "napasco_metric", "Watson-Crick", "provenance"):
            self.assertIn(present, self.raw, f"{present!r} expected in source mmCIF")
        # ... but the graph JSON keeps none of them (only identity + LW + backbone)
        for lost in ("napair_rmsd", "napasco", "nearest_curated", "Watson-Crick", "provenance", "orientation"):
            self.assertNotIn(lost, self.js, f"{lost!r} should be dropped from the JSON graph")

    def test_no_coordinates_when_source_has_none(self):
        self.assertNotIn("_atom_site", self.raw)          # annotation-only file: no coordinates
        for _, d in self.reloaded.nodes(data=True):
            self.assertNotIn("xyz_C1p", d)

    # --- 3. networkx build + visualization ----------------------------------
    def test_networkx_build_and_draw(self):
        G = self.reloaded
        S = nx.Graph()                                    # undirected view for drawing
        S.add_nodes_from(G.nodes(data=True))
        S.add_edges_from((u, v) for u, v in G.edges())
        self.assertEqual(S.number_of_nodes(), 26)

        pos = nx.spring_layout(S, seed=0)
        self.assertEqual(len(pos), S.number_of_nodes())

        out = Path(self.tmp.name) / "viz.png"
        fig, ax = plt.subplots()
        nx.draw(S, pos, ax=ax, node_size=60,
                labels={n: G.nodes[n]["nt_code"] for n in S}, font_size=6)
        fig.savefig(out)
        plt.close(fig)
        self.assertTrue(out.exists() and out.stat().st_size > 0)


@unittest.skipUnless(_find_dnatco_cif() is not None,
                     "no local DNATCO coordinate cif (run the dnatco fetch first)")
class TestDnatcoCoordsToJson(unittest.TestCase):
    """Coordinate mmCIF -> graph: representative atoms kept, full atom detail lost."""

    @classmethod
    def setUpClass(cls):
        from dnatco_to_graph import dnatco_cif_to_graph
        cls.cif = _find_dnatco_cif()
        cls.tmp = tempfile.TemporaryDirectory()
        cls.G = dnatco_cif_to_graph(cls.cif)
        p = Path(cls.tmp.name) / "g.json"
        dump_json(str(p), cls.G)
        cls.js = p.read_text()
        cls.reloaded = load_graph(str(p))

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_representative_coords_used(self):
        nodes = list(self.reloaded.nodes(data=True))
        with_c1 = [d for _, d in nodes if "xyz_C1p" in d]
        self.assertGreaterEqual(len(with_c1), 1)
        d = with_c1[0]
        self.assertEqual(len(d["xyz_C1p"]), 3)            # an (x, y, z) coordinate
        self.assertIn("xyz_glyN", d)                      # glycosidic N: present on every nucleotide
        # phosphate is absent on the 5'-terminal residue, so just require it somewhere
        self.assertTrue(any("xyz_P" in dd for _, dd in nodes))

    def test_full_atom_detail_lost(self):
        self.assertIn("xyz_C1p", self.js)                 # representative atom kept
        # the dozens of other atoms / the raw _atom_site table are not carried over
        for lost in ("_atom_site", "OP1", "O2'", "C2'"):
            self.assertNotIn(lost, self.js, f"{lost!r} should not be in the graph JSON")


if __name__ == "__main__":
    unittest.main(verbosity=2)
