#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import sys
import unittest

import awkward as ak
import numpy as np


MODULE_PATH = Path(__file__).with_name("plot_oo_v3_dijet_angles.py")
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("plot_oo_v3_dijet_angles", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
plotter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = plotter
SPEC.loader.exec_module(plotter)


class DijetAngleTest(unittest.TestCase):
    def test_extract_dijet_sorts_corrected_pt_and_wraps_delta_phi(self) -> None:
        arrays = ak.Array(
            {
                "jet4Pt": [[25.0, 40.0, 35.0], [31.0, 19.0], [45.0]],
                "jet4Eta": [[0.2, -1.0, 0.5], [0.1, -0.1], [0.0]],
                "jet4Phi": [[0.0, 3.0, -3.0], [0.2, 2.9], [1.0]],
            }
        )
        result = plotter.extract_dijet(arrays, 4, 30.0, 20.0)

        np.testing.assert_array_equal(result["selected"], [1.0, 0.0, 0.0])
        np.testing.assert_array_equal(result["event_indices"], [0])
        np.testing.assert_allclose(result["abs_delta_eta"], [1.5])
        np.testing.assert_allclose(result["delta_phi"], [2.0 * math.pi - 6.0])

    def test_weighted_quantiles_keep_fixed_physical_bounds(self) -> None:
        edges = plotter.weighted_quantile_edges(
            np.array([0.1, 0.2, 0.3, 0.4]),
            np.ones(4),
            lower=0.0,
            upper=1.0,
            bins=2,
        )
        self.assertEqual(edges[0], 0.0)
        self.assertEqual(edges[-1], 1.0)
        self.assertTrue(np.all(np.diff(edges) > 0.0))

    def test_normalized_histogram_has_unit_area(self) -> None:
        weights = np.array([1.0, 2.0, 3.0])
        selected = np.array([1.0, 1.0, 0.0])
        density, error, matrix, underflow, overflow = plotter.normalized_histogram(
            np.array([0.25, 1.25]),
            np.array([0, 1]),
            selected,
            weights,
            np.array([0.0, 0.5, 1.5]),
        )
        self.assertAlmostEqual(float(np.sum(density * np.array([0.5, 1.0]))), 1.0)
        self.assertTrue(np.all(np.isfinite(error)))
        self.assertEqual(matrix.shape, (3, 2))
        self.assertEqual((underflow, overflow), (0, 0))


if __name__ == "__main__":
    unittest.main()
