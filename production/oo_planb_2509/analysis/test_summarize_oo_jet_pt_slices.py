#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


MODULE_PATH = Path(__file__).with_name("summarize_oo_jet_pt_slices.py")
SPEC = importlib.util.spec_from_file_location("summarize_oo_jet_pt_slices", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
summary = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = summary
SPEC.loader.exec_module(summary)


def metadata(
    low: float,
    high: float | None,
    values: tuple[float, float, float, float, float, float, float, float],
):
    common = {
        "inputRoot": "/sample.root",
        "inputBytes": 123,
        "pairCount": 10,
        "uniqueSeedCount": 10,
        "weightSum": 8.0,
        "weightedSigmaSum": 16.0,
        "sigmaMergedMb": 2.0,
        "crossSectionFactor": 0.25,
        "effectiveEventCountFromWeights": 6.0,
        "normalization": "merged",
        "uncertainty": "paired jackknife",
        "ptMinGeV": low,
        "ptMaxGeV": high,
    }
    common["radii"] = {
        radius: {
            "variants": {
                "noPrehydro": {"integratedJetCrossSectionMb": values[index]},
                "withPrehydro": {"integratedJetCrossSectionMb": values[index + 1]},
            }
        }
        for index, radius in zip(range(0, len(values), 2), summary.RADII)
    }
    return common


class PtSliceSummaryTest(unittest.TestCase):
    def setUp(self) -> None:
        self.inclusive = metadata(
            30.0, None, (2.5, 2.25, 5.0, 4.5, 10.0, 9.0, 20.0, 18.0)
        )
        self.slices = [
            metadata(30.0, 50.0, (1.75, 1.5, 3.5, 3.0, 7.0, 6.0, 14.0, 12.0)),
            metadata(50.0, 80.0, (0.5, 0.5, 1.0, 1.0, 2.0, 2.0, 4.0, 4.0)),
            metadata(80.0, None, (0.25, 0.25, 0.5, 0.5, 1.0, 1.0, 2.0, 2.0)),
        ]

    def test_valid_partition_closes_to_inclusive(self) -> None:
        result = summary.validate_slices(self.slices, self.inclusive)
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(set(result["crossSectionClosureDifferenceMb"]), set(summary.RADII))
        self.assertEqual(summary.range_label(self.slices[-1]), ">80")

    def test_boundary_gap_is_rejected(self) -> None:
        self.slices[1]["ptMinGeV"] = 51.0
        with self.assertRaisesRegex(ValueError, "gap or overlap"):
            summary.validate_slices(self.slices, self.inclusive)

    def test_cross_section_mismatch_is_rejected(self) -> None:
        self.slices[0]["radii"]["0.4"]["variants"]["noPrehydro"][
            "integratedJetCrossSectionMb"
        ] = 7.1
        with self.assertRaisesRegex(ValueError, "do not close"):
            summary.validate_slices(self.slices, self.inclusive)


if __name__ == "__main__":
    unittest.main()
