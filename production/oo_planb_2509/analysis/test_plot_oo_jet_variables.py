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


if __name__ == "__main__":
    unittest.main()
