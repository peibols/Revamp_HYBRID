#!/usr/bin/env python3

from __future__ import annotations

import io
import math
import unittest

import analyze_oo_prehydro_pair as analysis


class ParsePythiaRunTest(unittest.TestCase):
    def test_reconstructs_weight_sum_and_uses_final_sigma_gen(self) -> None:
        stream = io.BytesIO(
            b"""# event 0
weight 2 cross 10 X 0 Y 0
3 4 0 0.13957 211 0
end
# event 1
weight 0.5 cross 12 X 0 Y 0
6 8 0 0.49368 321 0
end
"""
        )

        run = analysis.parse_pythia_run(stream, [0.0, 6.0, 12.0], 1.0)

        self.assertEqual(run.event_count, 2)
        self.assertEqual(run.weight_sum, 2.5)
        self.assertEqual(run.sigma_gen, 12.0)
        self.assertEqual(run.histogram, [2.0 / 6.0, 0.5 / 6.0])

    def test_applies_negative_wake_label(self) -> None:
        stream = io.BytesIO(
            b"""# event 0
weight 3 cross 8 X 0 Y 0
3 4 0 0.13957 211 0
3 4 0 0.13957 211 2
end
"""
        )

        run = analysis.parse_pythia_run(stream, [0.0, 10.0], 1.0)

        self.assertEqual(run.histogram, [0.0])

    def test_rejects_event_without_sigma_gen(self) -> None:
        stream = io.BytesIO(
            b"""# event 0
weight 1
3 4 0 0.13957 211 0
end
"""
        )

        with self.assertRaisesRegex(ValueError, "weight or sigmaGen"):
            analysis.parse_pythia_run(stream, [0.0, 10.0], 1.0)


class PythiaAggregateTest(unittest.TestCase):
    def test_matches_pythia_parallel_weighted_sigma_rule(self) -> None:
        aggregate = analysis.PythiaAggregate(1)
        aggregate.add(analysis.PythiaRun([2.0], 10.0, 2.0, 1))
        aggregate.add(analysis.PythiaRun([3.0], 20.0, 3.0, 1))

        stats = aggregate.stats()

        self.assertEqual(stats.event_count, 2)
        self.assertEqual(stats.run_count, 2)
        self.assertEqual(stats.weight_sum, 5.0)
        self.assertEqual(stats.sigma_gen, 16.0)
        self.assertEqual(stats.values, [16.0])
        self.assertTrue(math.isclose(stats.standard_errors[0], 5.0))


if __name__ == "__main__":
    unittest.main()
