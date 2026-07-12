#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

import awkward as ak
import numpy as np
import uproot


MODULE_PATH = Path(__file__).with_name("plot_oo_jet_paired_substructure.py")
SPEC = importlib.util.spec_from_file_location(
    "plot_oo_jet_paired_substructure", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
analysis = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = analysis
SPEC.loader.exec_module(analysis)


class PairedSubstructureTest(unittest.TestCase):
    def test_matched_flat_does_not_index_unselected_empty_event(self) -> None:
        values = ak.Array([[10.0], []])
        indices = ak.Array([[0], [0]])
        selection = ak.Array([[True], [False]])
        np.testing.assert_array_equal(
            analysis.matched_flat(values, indices, selection),
            np.array([10.0]),
        )

    def test_weighted_ratio_and_block_jackknife(self) -> None:
        estimate = analysis.mean_estimate(
            values=np.array([1.0, 3.0, 5.0, 7.0]),
            mask=np.ones(4, dtype=bool),
            event_indices=np.arange(4),
            event_weights=np.array([1.0, 2.0, 1.0, 2.0]),
            event_hydro_indices=np.array([0, 0, 1, 1]),
        )
        self.assertAlmostEqual(estimate.value, 26.0 / 6.0)
        self.assertEqual(estimate.raw_entries, 4)
        self.assertEqual(estimate.raw_events, 4)
        self.assertTrue(np.isfinite(estimate.event_error))
        self.assertTrue(np.isfinite(estimate.hydro_error))
        shifted = analysis.mean_estimate(
            values=np.array([2.0, 4.0, 6.0, 8.0]),
            mask=np.ones(4, dtype=bool),
            event_indices=np.arange(4),
            event_weights=np.array([1.0, 2.0, 1.0, 2.0]),
            event_hydro_indices=np.array([0, 0, 1, 1]),
        )
        difference = analysis.difference_estimate(shifted, estimate)
        self.assertAlmostEqual(difference.value, 1.0)
        self.assertAlmostEqual(difference.event_error, 0.0)
        self.assertAlmostEqual(difference.hydro_error, 0.0)

    def test_end_to_end_matched_selection_and_transitions(self) -> None:
        pair_ids = np.arange(8, dtype=np.uint64)
        weights = np.array([1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        no_pt = [[40.0], [42.0], [48.0], [52.0], [60.0], [65.0], [90.0], [100.0]]
        pre_pt = [[25.0], [39.0], [45.0], [49.0], [55.0], [59.0], [80.0], [88.0]]
        eta = [[0.1], [-0.2], [0.3], [-0.4], [0.2], [-0.1], [0.5], [-0.6]]
        no_sd = [[0], [0], [1], [1], [0], [1], [0], [1]]
        pre_sd = [[0], [1], [0], [1], [1], [1], [0], [0]]
        no_ptd = [[0.40], [0.42], [0.44], [0.46], [0.48], [0.50], [0.52], [0.54]]
        pre_ptd = [[value[0] + 0.01] for value in no_ptd]
        no_normal_ptd = [[value[0] + 0.05] for value in no_ptd]
        pre_normal_ptd = [[value[0] + 0.02] for value in no_normal_ptd]
        neff = [[1.4], [1.8], [4.0], [4.5], [9.0], [10.0], [14.0], [16.0]]
        hard_id = [[1], [2], [21], [21], [1], [2], [21], [21]]
        match_index = [[0] for _ in pair_ids]
        match_dr = [[0.02] for _ in pair_ids]

        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            root_path = work / "paired.root"
            no_branches: dict[str, object] = {
                "pairId": pair_ids,
                "eventWeight": weights,
            }
            pre_branches: dict[str, object] = {
                "pairId": pair_ids,
                "eventWeight": weights,
            }
            for radius_digit in analysis.RADIUS_DIGITS:
                prefix = f"jet{radius_digit}"
                no_branches[f"{prefix}Pt"] = ak.Array(no_pt)
                no_branches[f"{prefix}Eta"] = ak.Array(eta)
                no_branches[f"{prefix}PtD"] = ak.Array(no_ptd)
                no_branches[f"{prefix}NormalPtD"] = ak.Array(no_normal_ptd)
                no_branches[f"{prefix}NormalEffectiveMultiplicity"] = ak.Array(neff)
                no_branches[f"{prefix}SoftDropValid"] = ak.Array(no_sd)
                no_branches[f"{prefix}HardPartonId"] = ak.Array(hard_id)
                no_branches[f"{prefix}PairMatchIndex"] = ak.Array(match_index)
                no_branches[f"{prefix}PairMatchDR"] = ak.Array(match_dr)
                no_branches[f"{prefix}PairMatchOtherPt"] = ak.Array(pre_pt)
                no_branches[f"{prefix}PairMatchOtherHardPartonId"] = ak.Array(hard_id)

                pre_branches[f"{prefix}Pt"] = ak.Array(pre_pt)
                pre_branches[f"{prefix}PtD"] = ak.Array(pre_ptd)
                pre_branches[f"{prefix}NormalPtD"] = ak.Array(pre_normal_ptd)
                pre_branches[f"{prefix}SoftDropValid"] = ak.Array(pre_sd)
                pre_branches[f"{prefix}HardPartonId"] = ak.Array(hard_id)

            with uproot.recreate(root_path) as root_file:
                root_file["Pairs"] = {
                    "pairId": pair_ids,
                    "eventWeight": weights,
                    "seed": np.arange(100, 108, dtype=np.int64),
                    "hydroIndex": np.array(
                        [0, 0, 1, 1, 2, 2, 3, 3], dtype=np.int32
                    ),
                }
                root_file["noPrehydro/Jets"] = no_branches
                root_file["withPrehydro/Jets"] = pre_branches

            out_dir = work / "plots"
            metadata = analysis.run(
                argparse.Namespace(
                    input_root=root_path,
                    out_dir=out_dir,
                    pt_min=30.0,
                    pt_max=None,
                    abs_eta_max=2.0,
                    max_pair_match_dr_fraction=0.5,
                    prefix="test_paired",
                )
            )

            self.assertEqual(metadata["pairCount"], 8)
            self.assertEqual(metadata["uniqueHydroCount"], 4)
            self.assertEqual(metadata["radii"]["0.4"]["selectedMatchedJets"], 8)
            transition_path = out_dir / "test_paired_softdrop_transitions.tsv"
            with transition_path.open(newline="") as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            all_r04 = {
                row["transition"]: row
                for row in rows
                if row["radius"] == "0.4" and row["flavor"] == "all"
            }
            expected = {
                "fail_to_fail": 2.0 / 9.0,
                "fail_to_pass": 3.0 / 9.0,
                "pass_to_fail": 2.0 / 9.0,
                "pass_to_pass": 2.0 / 9.0,
            }
            for transition, fraction in expected.items():
                self.assertAlmostEqual(
                    float(all_r04[transition]["fraction"]), fraction
                )
                self.assertEqual(
                    int(all_r04[transition]["rawTransitionJets"]), 2
                )
            self.assertAlmostEqual(
                float(all_r04["fail_to_pass"]["weightedTransitionJets"]),
                3.0,
            )

            fail_metadata = metadata["radii"]["0.4"]["flavors"]["all"][
                "softDropFailure"
            ]
            self.assertAlmostEqual(
                fail_metadata["no_fail_fraction"]["value"], 5.0 / 9.0
            )
            self.assertAlmostEqual(
                fail_metadata["pre_fail_fraction"]["value"], 4.0 / 9.0
            )
            self.assertAlmostEqual(
                fail_metadata["delta_fail_fraction"]["value"], -1.0 / 9.0
            )
            paired = metadata["radii"]["0.4"]["flavors"]["all"][
                "observables"
            ]
            self.assertAlmostEqual(
                paired["all"]["delta_ptd"]["value"], 0.01
            )
            self.assertAlmostEqual(
                paired["all"]["delta_normal_ptd"]["value"], 0.02
            )
            for radius_digit in analysis.RADIUS_DIGITS:
                self.assertTrue(
                    (
                        out_dir
                        / f"test_paired_R0{radius_digit}_paired_substructure.pdf"
                    ).is_file()
                )
            self.assertTrue(
                (out_dir / "test_paired_paired_substructure_summary.pdf").is_file()
            )
            self.assertTrue(
                (out_dir / "test_paired_paired_observables.tsv").is_file()
            )
            self.assertTrue(
                (out_dir / "test_paired_single_many_contrasts.tsv").is_file()
            )
            self.assertTrue((out_dir / "test_paired_metadata.json").is_file())


if __name__ == "__main__":
    unittest.main()
