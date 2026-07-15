#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest

import awkward as ak
import numpy as np


DIRECTORY = Path(__file__).parent


def load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, DIRECTORY / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


load("plot_oo_jet_variables", "plot_oo_jet_variables.py")
plotter = load("plot_oo_jet_formation_time", "plot_oo_jet_formation_time.py")


class FormationTimePlotTest(unittest.TestCase):
    def test_root_times_are_scaled_to_e_over_q2(self) -> None:
        self.assertEqual(plotter.stored_tau_f_scale("oo-paired-root-v5"), 1.0)
        self.assertEqual(plotter.stored_tau_f_scale("oo-paired-root-v6"), 1.0)
        self.assertEqual(plotter.stored_tau_f_scale("oo-paired-root-v7"), 0.5)
        self.assertEqual(plotter.stored_tau_f_scale("oo-paired-root-v8"), 1.0)
        with self.assertRaises(ValueError):
            plotter.stored_tau_f_scale("oo-paired-root-v4")

    def test_weighted_yield_applies_event_weight_once(self) -> None:
        weights = np.array([1.0, 2.0, 3.0])
        matrix = np.array([[1.0, 0.0], [2.0, 0.0], [0.0, 6.0]])
        edges = np.array([-1.0, 0.0, 2.0])
        values, errors = plotter.weighted_yield_from_matrix(matrix, weights, edges)

        np.testing.assert_allclose(values, [3.0 / 6.0, 6.0 / 6.0 / 2.0])
        self.assertTrue(np.all(np.isfinite(errors)))

    def test_selected_jet_normalization_closes_to_splits_per_jet(self) -> None:
        matrix = np.array([[1.0, 0.0], [0.0, 2.0], [3.0, 0.0]])
        weights = np.array([1.0, 2.0, 3.0])
        selected_counts = np.array([1.0, 1.0, 2.0])
        edges = np.array([0.0, 1.0, 2.0])
        values, errors = plotter.selected_jet_normalized_yield(
            matrix, selected_counts, weights, edges
        )
        expected_splits_per_jet = float(np.sum(matrix)) / float(
            np.sum(weights * selected_counts)
        )
        self.assertAlmostEqual(float(np.sum(values * np.diff(edges))), expected_splits_per_jet)
        self.assertTrue(np.all(np.isfinite(errors)))

    def test_corrected_pt_intervals_are_disjoint(self) -> None:
        arrays = ak.Array(
            {
                "radius": [0.4] * 7,
                "correctedPt": [30.0, 40.0, 50.0, 60.0, 80.0, 90.0, np.nan],
            }
        )
        low = plotter.jet_selection(arrays, 0.4, 30.0, 50.0)
        middle = plotter.jet_selection(arrays, 0.4, 50.0, 80.0)
        high = plotter.jet_selection(arrays, 0.4, 80.0, None)

        np.testing.assert_array_equal(
            low, [False, True, True, False, False, False, False]
        )
        np.testing.assert_array_equal(
            middle, [False, False, False, True, True, False, False]
        )
        np.testing.assert_array_equal(
            high, [False, False, False, False, False, True, False]
        )
        self.assertFalse(np.any((low & middle) | (low & high) | (middle & high)))

    def test_weighted_quantile_uses_split_event_weights(self) -> None:
        values = np.array([1.0, 2.0, 10.0])
        weights = np.array([1.0, 1.0, 8.0])
        median = plotter.weighted_quantile(values, weights, [0.5])[0]
        self.assertGreater(median, 5.0)

    def test_paired_ratio_of_ratios_uses_shared_runs(self) -> None:
        weights = np.array([1.0, 2.0, 3.0])
        value, error = plotter.paired_ratio_of_ratios(
            np.array([2.0, 1.0, 3.0]),
            np.array([1.0, 1.0, 2.0]),
            np.array([1.0, 2.0, 2.0]),
            np.array([1.0, 1.0, 1.0]),
            weights,
        )
        expected = ((2.0 + 2.0 + 9.0) / (1.0 + 2.0 + 6.0)) / (
            (1.0 + 4.0 + 6.0) / (1.0 + 2.0 + 3.0)
        )
        self.assertAlmostEqual(value, expected)
        self.assertTrue(np.isfinite(error))
        self.assertGreater(error, 0.0)


if __name__ == "__main__":
    unittest.main()
