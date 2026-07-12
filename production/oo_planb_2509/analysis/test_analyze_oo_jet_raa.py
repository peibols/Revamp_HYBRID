#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np


MODULE_PATH = Path(__file__).with_name("analyze_oo_jet_raa.py")
SPEC = importlib.util.spec_from_file_location("analyze_oo_jet_raa", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
analyzer = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = analyzer
SPEC.loader.exec_module(analyzer)


class JetRaaNormalizationTest(unittest.TestCase):
    def test_pythia_parallel_aggregate_and_jackknife(self) -> None:
        matrix = np.array([[1.0, 0.0], [0.0, 2.0], [3.0, 1.0]])
        weight_sums = np.array([1.0, 2.0, 3.0])
        sigma_gen = np.array([10.0, 20.0, 30.0])
        values, errors, weight_sum, sigma_merged = analyzer.aggregate_spectrum(
            matrix, weight_sums, sigma_gen
        )
        self.assertEqual(weight_sum, 6.0)
        self.assertAlmostEqual(sigma_merged, 140.0 / 6.0)
        np.testing.assert_allclose(values, (140.0 / 36.0) * np.array([4.0, 3.0]))
        self.assertTrue(np.all(np.isfinite(errors)))

    def test_independent_ratio_error(self) -> None:
        ratio, error = analyzer.independent_ratio(
            np.array([8.0]), np.array([0.8]), np.array([10.0]), np.array([0.5])
        )
        self.assertAlmostEqual(ratio[0], 0.8)
        self.assertAlmostEqual(error[0], 0.8 * math.sqrt(0.1**2 + 0.05**2))

    def test_streaming_helper_clusters_all_three_radii_and_bin_edges(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_text:
            temporary = Path(temporary_text)
            try:
                helper = analyzer.build_helper(
                    Path(__file__).with_name("oo_pp_jet_spectrum.cc"), temporary, "g++", True
                )
            except (FileNotFoundError, subprocess.CalledProcessError) as error:
                self.skipTest(f"FastJet compiler unavailable: {error}")
            output = temporary / "spectrum.tsv"
            stream = """# OORUN 7
# event 0
weight 2 cross 10 X 0 Y 0
40 0 0 0 211 0
-20 0 0 0 -211 0
100 0 0 0 21 -2
end
# OOENDRUN
"""
            subprocess.run(
                [
                    str(helper),
                    "--output",
                    str(output),
                    "--bins",
                    "20,40",
                    "--jet-abs-eta-max",
                    "2",
                ],
                input=stream,
                text=True,
                check=True,
            )
            with output.open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 3)
            for row in rows:
                self.assertEqual(int(row["run_id"]), 7)
                self.assertEqual(int(row["event_count"]), 1)
                self.assertEqual(int(row["raw_jet_count"]), 1)
                self.assertAlmostEqual(float(row["weighted_density"]), 0.1)


if __name__ == "__main__":
    unittest.main()
