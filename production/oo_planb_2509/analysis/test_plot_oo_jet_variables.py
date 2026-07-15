#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import sys
import unittest

import awkward as ak
import numpy as np


MODULE_PATH = Path(__file__).with_name("plot_oo_jet_variables.py")
SPEC = importlib.util.spec_from_file_location("plot_oo_jet_variables", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
plotter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = plotter
SPEC.loader.exec_module(plotter)


class JetWeightingTest(unittest.TestCase):
    def setUp(self) -> None:
        self.weights = np.array([1.0, 2.0, 3.0])
        self.sigma_gen = np.array([10.0, 20.0, 30.0])

    def test_pythia_parallel_factor(self) -> None:
        factor, weight_sum, sigma_merged = plotter.pythia_parallel_factor(
            self.weights, self.sigma_gen
        )
        self.assertEqual(weight_sum, 6.0)
        self.assertAlmostEqual(sigma_merged, 140.0 / 6.0)
        self.assertAlmostEqual(factor, 140.0 / 36.0)

    def test_differential_histogram_uses_one_merged_normalization(self) -> None:
        values = ak.Array([[0.2, 0.8], [0.4], []])
        result = plotter.differential_histogram(
            values,
            self.weights,
            self.sigma_gen,
            np.array([0.0, 0.5, 1.0]),
        )
        factor = 140.0 / 36.0
        np.testing.assert_allclose(result.values, [factor * 3.0 / 0.5, factor * 1.0 / 0.5])
        self.assertEqual(result.raw_entries, 3)
        self.assertEqual(result.weighted_entries, 4.0)
        self.assertTrue(np.all(np.isfinite(result.errors)))

    def test_integrated_cross_section_and_paired_ratio(self) -> None:
        no_counts = np.array([1.0, 2.0, 0.0])
        pre_counts = np.array([1.0, 1.0, 1.0])
        result = plotter.integrated_cross_section(no_counts, self.weights, self.sigma_gen)
        self.assertAlmostEqual(result.cross_section, (140.0 / 36.0) * 5.0)
        ratio, error = plotter.paired_integrated_ratio(pre_counts, no_counts, self.weights)
        self.assertAlmostEqual(ratio, 6.0 / 5.0)
        self.assertTrue(math.isfinite(error))
        self.assertGreater(error, 0.0)

    def test_histogram_ratio_uses_paired_weighted_contributions(self) -> None:
        denominator = np.array([[1.0, 0.0], [2.0, 1.0], [1.0, 3.0]])
        numerator = np.array([[1.0, 1.0], [1.0, 1.0], [2.0, 2.0]])
        ratio, error = plotter.paired_ratio(numerator, denominator)
        np.testing.assert_allclose(ratio, [1.0, 1.0])
        self.assertTrue(np.all(np.isfinite(error)))

    def test_shape_normalization_and_ratio_use_selected_jet_area(self) -> None:
        edges = np.array([0.0, 1.0, 2.0])
        no = plotter.differential_histogram(
            ak.Array([[0.5], [0.5, 1.5], [1.5]]),
            self.weights,
            self.sigma_gen,
            edges,
        )
        pre = plotter.differential_histogram(
            ak.Array([[0.5], [1.5], [1.5, 1.5]]),
            self.weights,
            self.sigma_gen,
            edges,
        )
        no_counts = np.array([1.0, 2.0, 1.0])
        pre_counts = np.array([1.0, 1.0, 2.0])
        plotter.normalize_histogram_to_selected_jets(
            no, no_counts, self.weights, edges
        )
        plotter.normalize_histogram_to_selected_jets(
            pre, pre_counts, self.weights, edges
        )
        self.assertAlmostEqual(float(np.sum(no.normalized_values)), 1.0)
        self.assertAlmostEqual(float(np.sum(pre.normalized_values)), 1.0)
        ratio, error = plotter.paired_normalized_ratio(
            pre.weighted_matrix,
            no.weighted_matrix,
            pre_counts,
            no_counts,
            self.weights,
        )
        np.testing.assert_allclose(ratio, pre.normalized_values / no.normalized_values)
        self.assertTrue(np.all(np.isfinite(error)))

    def test_pt_slices_partition_the_inclusive_selection(self) -> None:
        pt = np.array([30.0, 40.0, 50.0, 60.0, 80.0, 90.0, np.nan])
        low = plotter.jet_pt_selection(pt, 30.0, 50.0)
        middle = plotter.jet_pt_selection(pt, 50.0, 80.0)
        high = plotter.jet_pt_selection(pt, 80.0, None)
        inclusive = plotter.jet_pt_selection(pt, 30.0, None)

        np.testing.assert_array_equal(low, [False, True, True, False, False, False, False])
        np.testing.assert_array_equal(middle, [False, False, False, True, True, False, False])
        np.testing.assert_array_equal(high, [False, False, False, False, False, True, False])
        np.testing.assert_array_equal(low | middle | high, inclusive)
        self.assertFalse(np.any((low & middle) | (low & high) | (middle & high)))

    def test_bounded_pt_histogram_uses_the_requested_edges(self) -> None:
        spec = plotter.variable_specs(4, 30.0, 50.0)[0]
        np.testing.assert_allclose(spec.edges, np.linspace(30.0, 50.0, 11))
        self.assertEqual(spec.xscale, "linear")
        self.assertEqual(plotter.pt_range_label(30.0, 50.0), r"$30<p_T^{\rm jet}\leq 50$ GeV")

    def test_radius_01_and_02_have_dedicated_shape_ranges(self) -> None:
        radius_01 = {spec.key: spec for spec in plotter.variable_specs(1, 20.0)}
        self.assertAlmostEqual(float(radius_01["rg"].edges[-1]), 0.125)
        self.assertAlmostEqual(float(radius_01["girth"].edges[-1]), 0.08)
        specs = {spec.key: spec for spec in plotter.variable_specs(2, 20.0)}
        self.assertAlmostEqual(float(specs["rg"].edges[-1]), 0.25)
        self.assertAlmostEqual(float(specs["girth"].edges[-1]), 0.16)
        self.assertEqual(plotter.RADIUS_DIGITS, (1, 2, 4, 8))

    def test_signed_total_multiplicity_has_negative_bins(self) -> None:
        for radius in plotter.RADIUS_DIGITS:
            specs = {spec.key: spec for spec in plotter.variable_specs(radius, 20.0)}
            self.assertIn("totalmult", specs)
            self.assertEqual(specs["totalmult"].branch_suffix, "TotalMult")
            self.assertLess(float(specs["totalmult"].edges[0]), 0.0)
            self.assertGreater(float(specs["totalmult"].edges[-1]), 0.0)

    def test_factor_two_rebin_preserves_soft_drop_failure_bin(self) -> None:
        original = {
            spec.key: spec for spec in plotter.variable_specs(4, 30.0)
        }
        rebinned = {
            spec.key: spec
            for spec in plotter.variable_specs(
                4, 30.0, substructure_rebin_factor=2
            )
        }
        for key in ("zg", "rg"):
            self.assertEqual(rebinned[key].edges[0], original[key].edges[0])
            self.assertEqual(rebinned[key].edges[1], original[key].edges[1])
            self.assertEqual(
                len(rebinned[key].edges) - 2,
                (len(original[key].edges) - 2) // 2,
            )
        self.assertEqual(len(rebinned["ptd"].edges) - 1, 10)
        np.testing.assert_array_equal(rebinned["pt"].edges, original["pt"].edges)

    def test_factor_two_rebin_retains_odd_trailing_bin(self) -> None:
        original = {
            spec.key: spec for spec in plotter.variable_specs(4, 30.0)
        }["maxkt"].edges
        rebinned = {
            spec.key: spec
            for spec in plotter.variable_specs(
                4, 30.0, substructure_rebin_factor=2
            )
        }["maxkt"].edges
        self.assertEqual(len(original) - 1, 19)
        self.assertEqual(len(rebinned) - 1, 10)
        self.assertEqual(rebinned[-2], original[-2])
        self.assertEqual(rebinned[-1], original[-1])

    def test_soft_drop_specs_reserve_equal_width_failure_bin(self) -> None:
        specs = {spec.key: spec for spec in plotter.variable_specs(4, 30.0)}
        for key, physical_minimum in (("zg", 0.1), ("rg", 0.0)):
            spec = specs[key]
            self.assertTrue(spec.soft_drop_failure_bin)
            self.assertAlmostEqual(float(spec.edges[1]), physical_minimum)
            self.assertAlmostEqual(
                float(spec.edges[1] - spec.edges[0]),
                float(spec.edges[2] - spec.edges[1]),
            )

    def test_failed_soft_drop_jets_fill_first_bin(self) -> None:
        spec = {spec.key: spec for spec in plotter.variable_specs(4, 30.0)}["zg"]
        values = ak.Array([[math.nan, 0.2], [math.nan, 0.3], []])
        valid = ak.Array([[False, True], [False, True], []])
        encoded = plotter.soft_drop_histogram_values(values, valid, spec)
        result = plotter.differential_histogram(
            encoded, self.weights, self.sigma_gen, spec.edges
        )

        factor = 140.0 / 36.0
        first_bin_width = float(spec.edges[1] - spec.edges[0])
        self.assertEqual(result.raw_entries, 4)
        self.assertEqual(result.nonfinite_entries, 0)
        self.assertAlmostEqual(result.values[0] * first_bin_width, factor * 3.0)
        self.assertAlmostEqual(
            float(np.sum(result.values * np.diff(spec.edges))), factor * 6.0
        )

        valid_only_mean = plotter.weighted_mean(values[valid], self.weights)
        self.assertEqual(valid_only_mean.raw_entries, 2)
        self.assertAlmostEqual(valid_only_mean.mean, 0.8 / 3.0)


if __name__ == "__main__":
    unittest.main()
