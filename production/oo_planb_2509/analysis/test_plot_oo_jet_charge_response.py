#!/usr/bin/env python3

from __future__ import annotations

import argparse
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import awkward as ak
import numpy as np
import uproot


MODULE_PATH = Path(__file__).with_name("plot_oo_jet_charge_response.py")
SPEC = importlib.util.spec_from_file_location("plot_oo_jet_charge_response", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
analysis = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = analysis
SPEC.loader.exec_module(analysis)


class ChargeResponseTest(unittest.TestCase):
    def test_flavor_masks(self) -> None:
        ids = np.array([1, -2, 6, 21, 0, 22])
        np.testing.assert_array_equal(
            analysis.flavor_mask(ids, "quark"),
            np.array([True, True, True, False, False, False]),
        )
        np.testing.assert_array_equal(
            analysis.flavor_mask(ids, "gluon"),
            np.array([False, False, False, True, False, False]),
        )
        np.testing.assert_array_equal(
            analysis.flavor_mask(ids, "quark", np.array([1, 2, 6, 21, 0, 22])),
            np.array([True, False, True, False, False, False]),
        )

    def test_binned_fractional_response(self) -> None:
        result = analysis.binned_response(
            proxy=np.array([1.0, 1.5, 4.0, 4.5]),
            fractional_shift=np.array([0.1, 0.2, 0.3, 0.4]),
            baseline_pt=np.array([10.0, 10.0, 20.0, 20.0]),
            event_indices=np.array([0, 1, 2, 3]),
            event_weights=np.ones(4),
            edges=np.array([0.0, 2.0, 5.0]),
        )
        np.testing.assert_allclose(result.mean, np.array([0.15, 0.35]))
        np.testing.assert_allclose(
            result.pt_weighted_fraction, np.array([0.15, 0.35])
        )
        np.testing.assert_array_equal(result.raw_entries, np.array([2, 2]))
        self.assertTrue(np.all(np.isfinite(result.mean_error)))

    def test_end_to_end_outputs(self) -> None:
        pair_ids = np.arange(8, dtype=np.uint64)
        weights = np.array([1.0, 1.4, 0.8, 1.2, 1.1, 0.9, 1.3, 0.7])
        pt = [[40.0], [42.0], [48.0], [52.0], [60.0], [65.0], [90.0], [100.0]]
        eta = [[0.1], [-0.2], [0.3], [-0.4], [0.2], [-0.1], [0.5], [-0.6]]
        neff = [[1.4], [1.8], [4.0], [4.5], [9.0], [10.0], [14.0], [16.0]]
        leading = [[0.9], [0.85], [0.6], [0.55], [0.35], [0.3], [0.22], [0.2]]
        nsd = [[0], [0], [1], [1], [2], [2], [3], [3]]
        max_kt = [[float("nan")], [0.7], [1.5], [2.5], [5.0], [7.0], [12.0], [18.0]]
        hard_id = [[1], [2], [21], [21], [1], [2], [21], [21]]
        match_index = [[0] for _ in pair_ids]
        match_dr = [[0.02] for _ in pair_ids]
        other_pt = [[38.0], [40.0], [46.0], [49.0], [54.0], [58.0], [80.0], [88.0]]

        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            root_path = work / "paired.root"
            branches: dict[str, object] = {
                "pairId": pair_ids,
                "eventWeight": weights,
            }
            for radius_digit in analysis.RADIUS_DIGITS:
                prefix = f"jet{radius_digit}"
                branches[f"{prefix}Pt"] = ak.Array(pt)
                branches[f"{prefix}Eta"] = ak.Array(eta)
                branches[f"{prefix}NormalEffectiveMultiplicity"] = ak.Array(neff)
                branches[f"{prefix}LeadingNormalFraction"] = ak.Array(leading)
                branches[f"{prefix}NSD"] = ak.Array(nsd)
                branches[f"{prefix}MaxKt"] = ak.Array(max_kt)
                branches[f"{prefix}HardPartonId"] = ak.Array(hard_id)
                branches[f"{prefix}PairMatchIndex"] = ak.Array(match_index)
                branches[f"{prefix}PairMatchDR"] = ak.Array(match_dr)
                branches[f"{prefix}PairMatchOtherPt"] = ak.Array(other_pt)
                branches[f"{prefix}PairMatchOtherHardPartonId"] = ak.Array(hard_id)
            with uproot.recreate(root_path) as root_file:
                root_file["Pairs"] = {
                    "pairId": pair_ids,
                    "eventWeight": weights,
                    "seed": np.arange(100, 108, dtype=np.int64),
                }
                root_file["noPrehydro/Jets"] = branches

            out_dir = work / "plots"
            metadata = analysis.run(
                argparse.Namespace(
                    input_root=root_path,
                    out_dir=out_dir,
                    pt_min=30.0,
                    pt_max=None,
                    abs_eta_max=2.0,
                    max_pair_match_dr_fraction=0.5,
                    prefix="test_charge",
                )
            )
            self.assertEqual(metadata["pairCount"], 8)
            for radius_digit in analysis.RADIUS_DIGITS:
                self.assertTrue(
                    (out_dir / f"test_charge_R0{radius_digit}_charge_response.pdf").is_file()
                )
                self.assertTrue(
                    (out_dir / f"test_charge_R0{radius_digit}_charge_response.png").is_file()
                )
            self.assertTrue((out_dir / "test_charge_charge_response.tsv").is_file())
            self.assertTrue((out_dir / "test_charge_single_many_summary.tsv").is_file())
            self.assertTrue(
                (out_dir / "test_charge_quark_single_many_summary.pdf").is_file()
            )
            self.assertTrue(
                (out_dir / "test_charge_quark_single_many_summary.png").is_file()
            )
            self.assertTrue(
                (out_dir / "test_charge_charge_response_metadata.json").is_file()
            )


if __name__ == "__main__":
    unittest.main()
