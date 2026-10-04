#!/usr/bin/env python3
"""Checks selection, exact accounting, risk placement and identical HTML pixels."""
import base64
import copy
import json
import math
import re
import tempfile
import unittest
from pathlib import Path
from dataclasses import replace
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from compass_os.sample_path_render import (          # noqa: E402
    _demo, Node, Style, select_nodes, plot_sample_path, aggregate_cox_predictor)


class RendererChecks(unittest.TestCase):
    def test_five_model_budgets(self):
        for model, multi, concepts_expected in [("M1", False, 16), ("M2", False, 13),
                ("M2", True, 12), ("M3", False, 12), ("M3", True, 11)]:
            d = _demo(model)
            if multi:
                d.cancer_types = ["A", "B"]
                d.predictors = list(d.predictors) + [Node("Cancer type", .1, "cancer", .1)]
                d.decomposition = []
            selected = select_nodes(d)
            self.assertEqual(selected["counts"]["predictors"], 16)
            self.assertEqual(sum(n.kind == "concept" for n in selected["predictors"]), concepts_expected)
            self.assertEqual(sum(n.kind == "pc" for n in selected["predictors"]), model == "M3")

    def test_threshold_and_matching_gene_columns(self):
        d = _demo("M3")
        selected = select_nodes(d, Style(gene_threshold=.85, max_genes=8))
        self.assertLessEqual(len(selected["gene_score"]), 8)
        self.assertTrue(all(abs(n.score) >= .85 for n in selected["gene_score"]))
        self.assertEqual([n.name for n in selected["gene_tpm"]], [n.name for n in selected["gene_score"]])

    def test_bookkeeping_removed(self):
        d = _demo("M3")
        d.predictors = list(d.predictors) + [Node("Other concepts", 999, "concept", 999),
                                             Node("Model standardization offset", 999, "concept", 999)]
        d.decomposition = []
        self.assertFalse(any(n.score == 999 for n in select_nodes(d)["predictors"]))

    def test_empty_threshold_is_honest(self):
        selected = select_nodes(_demo("M3"), Style(gene_threshold=100))
        self.assertEqual(selected["gene_score"], [])
        self.assertEqual(selected["gene_tpm"], [])
        self.assertEqual(selected["gene_edges"], [])

    def test_missing_tpm_fails(self):
        d = _demo("M2"); d.gene_tpm = []
        with self.assertRaisesRegex(ValueError, "Missing expression"):
            select_nodes(d)

    def test_missing_pc_fails(self):
        d = _demo("M3"); d.predictors = [n for n in d.predictors if n.kind != "pc"]; d.decomposition = []
        with self.assertRaisesRegex(ValueError, "aggregated predictors"):
            select_nodes(d)

    def test_signed_decomposition_must_equal_eta(self):
        d = _demo("M3"); d.risk += 1
        with self.assertRaisesRegex(ValueError, "not eta"):
            select_nodes(d)

    def test_encoded_column_sum(self):
        node = aggregate_cox_predictor("PC1–PC10", ["PC1", "PC2"], {"PC1": .4, "PC2": -.3}, kind="pc")
        self.assertTrue(math.isclose(node.contribution, .1))
        self.assertEqual(node.score, node.contribution)
        with self.assertRaises(ValueError):
            aggregate_cox_predictor("Stage", ["Stage", "Stage"], {"Stage": 1}, kind="clinical")

    def test_risk_endpoints_and_nonmedian_cutoff(self):
        d = _demo("M2"); d.decomposition = []
        d.cohort_risks = [0, 1, 2, 3, 4]; d.rank = 1; d.risk = 0; d.percentile = 0; d.cutoff = 3
        r = plot_sample_path(d)
        self.assertEqual(r.cutoff_percentile, 70)
        from matplotlib.patches import Polygon
        marker = next(p for p in r.figure.axes[0].patches if isinstance(p, Polygon))
        self.assertAlmostEqual(marker.get_xy()[2][1], 17)
        r.close()
        d.percentile = 100; d.rank = 5; d.risk = 4
        r = plot_sample_path(d)
        marker = next(p for p in r.figure.axes[0].patches if isinstance(p, Polygon))
        self.assertAlmostEqual(marker.get_xy()[2][1], 80)
        r.close()

    def test_median_with_ties(self):
        d = _demo("M2"); d.cohort_risks = [1, 1, 1]; d.rank = 2; d.decomposition = []
        r = plot_sample_path(d)
        self.assertEqual(r.cutoff, 1)
        self.assertEqual(r.cutoff_percentile, 50)
        r.close()

    def test_html_is_same_png_and_node_values_unchanged(self):
        d = _demo("M3"); original = copy.deepcopy(d)
        r = plot_sample_path(d)
        with tempfile.TemporaryDirectory() as td:
            files = r.save(Path(td)/"M3")
            embedded = re.search(r"base64,([^\"]+)", files["html"].read_text()).group(1)
            self.assertEqual(base64.b64decode(embedded), files["png"].read_bytes())
            self.assertTrue(files["pdf"].read_bytes().startswith(b"%PDF"))
            info = json.loads(files["selection.json"].read_text())
            self.assertEqual(info["percentile"], d.percentile)
            self.assertEqual(info["risk"], d.risk)
        self.assertEqual(d, original)
        r.close()

    def test_single_style_for_all_connections(self):
        r = plot_sample_path(_demo("M3"))
        from matplotlib.patches import PathPatch
        paths = [p for p in r.figure.axes[0].patches if isinstance(p, PathPatch)]
        self.assertTrue(paths)
        self.assertEqual(len({p.get_linewidth() for p in paths}), 1)
        self.assertEqual(len({p.get_edgecolor() for p in paths}), 1)
        r.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
