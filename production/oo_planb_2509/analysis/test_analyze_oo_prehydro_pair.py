#!/usr/bin/env python3

from __future__ import annotations

import csv
import io
import math
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

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


class PlotBinningTest(unittest.TestCase):
    def test_uses_fixed_equal_precision_bins(self) -> None:
        self.assertEqual(
            analysis.DEFAULT_BINS,
            [4, 5, 7, 10, 14, 24, 36, 50, 80, 150],
        )

    def test_uses_geometric_center_for_log_axis(self) -> None:
        self.assertEqual(analysis.logarithmic_bin_center(4.0, 9.0), 6.0)
        with self.assertRaises(ValueError):
            analysis.logarithmic_bin_center(0.0, 1.0)


class ShiftedManifestTest(unittest.TestCase):
    @staticmethod
    def write_manifest(path: Path, task_ids: list[int]) -> None:
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "task_id",
                    "hard_seed",
                    "hydro_slot",
                    "hydro_event_id",
                    "hydro_ncoll",
                    "hydro_payload_sha256",
                ],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for task_id in task_ids:
                writer.writerow(
                    {
                        "task_id": task_id,
                        "hard_seed": 900_000 + task_id,
                        "hydro_slot": 0,
                        "hydro_event_id": 1,
                        "hydro_ncoll": 2,
                        "hydro_payload_sha256": "a" * 64,
                    }
                )

    def test_accepts_shifted_contiguous_task_range(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.tsv"
            self.write_manifest(path, [50_000, 50_001])
            assignments = analysis.load_aa_task_manifest(path)
            self.assertEqual(set(assignments), {50_000, 50_001})

    def test_rejects_shifted_manifest_gap(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.tsv"
            self.write_manifest(path, [50_000, 50_002])
            with self.assertRaisesRegex(ValueError, "one contiguous range"):
                analysis.load_aa_task_manifest(path)


class MultipleAaSourceTest(unittest.TestCase):
    event = b"""# event 0
weight 1 cross 2
4.5 0 0 0.13957 211 0
end
"""

    @staticmethod
    def add_member(archive: tarfile.TarFile, name: str, payload: bytes) -> None:
        member = tarfile.TarInfo(name)
        member.size = len(payload)
        archive.addfile(member, io.BytesIO(payload))

    def make_snapshot(self, root: Path, kind: str) -> None:
        status_dir = root / "status" / kind
        output_dir = root / "outputs" / kind
        status_dir.mkdir(parents=True)
        output_dir.mkdir(parents=True)
        (status_dir / "chunk_0.txt").write_text("status=success\n")
        with tarfile.open(output_dir / "chunk_0.tar.gz", "w:gz") as archive:
            if kind == "aa":
                self.add_member(
                    archive, "runs/task_00000/HYBRID_Hadrons.out", self.event
                )
                self.add_member(
                    archive,
                    "runs/task_00000_prehydro/HYBRID_Hadrons.out",
                    self.event,
                )
            else:
                self.add_member(archive, "runs/pp/HYBRID_Hadrons.out", self.event)

    def test_accepts_distinct_sources_and_rejects_duplicates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            primary = Path(tmp) / "primary"
            continuation = Path(tmp) / "continuation"
            primary.mkdir()
            continuation.mkdir()
            self.assertEqual(
                analysis.distinct_aa_sources(primary, [continuation]),
                [primary, continuation],
            )
            with self.assertRaisesRegex(ValueError, "must be distinct"):
                analysis.distinct_aa_sources(primary, [primary])

    def test_merges_overlapping_chunk_ids_from_distinct_sources(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            primary = root / "aa10k"
            continuation = root / "aa20k"
            pp = root / "pp1m"
            output = root / "analysis"
            self.make_snapshot(primary, "aa")
            self.make_snapshot(continuation, "aa")
            self.make_snapshot(pp, "pp")

            argv = [
                "analyzer",
                "--local-eos", str(primary),
                "--additional-aa-local-eos", str(continuation),
                "--pp-local-eos", str(pp),
                "--out-dir", str(output),
                "--require-paired-aa",
                "--aa-events-per-chunk", "1",
            ]
            with patch.object(sys, "argv", argv), patch.object(
                analysis, "maybe_plot"
            ):
                self.assertEqual(analysis.main(), 0)

            with (output / "aa_sources.tsv").open() as handle:
                source_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["chunks_used"] for row in source_rows], ["1", "1"])
            with (output / "oo5360_c0_5_prehydro_overlay_raa.tsv").open() as handle:
                raa_rows = list(csv.DictReader(handle, delimiter="\t"))
            first_bin = {
                row["variant"]: row
                for row in raa_rows
                if row["pt_low"] == "4" and row["pt_high"] == "5"
            }
            self.assertEqual(set(first_bin), set(analysis.AA_VARIANTS))
            self.assertTrue(all(row["aa_events"] == "2" for row in first_bin.values()))
            self.assertTrue(all(row["aa_runs"] == "2" for row in first_bin.values()))

    def test_skips_success_archive_with_no_complete_event(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            valid = root / "valid"
            invalid = root / "invalid"
            pp = root / "pp"
            output = root / "analysis"
            self.make_snapshot(valid, "aa")
            self.make_snapshot(invalid, "aa")
            self.make_snapshot(pp, "pp")
            incomplete = b""
            with tarfile.open(
                invalid / "outputs/aa/chunk_0.tar.gz", "w:gz"
            ) as archive:
                self.add_member(
                    archive, "runs/task_00000/HYBRID_Hadrons.out", incomplete
                )
                self.add_member(
                    archive,
                    "runs/task_00000_prehydro/HYBRID_Hadrons.out",
                    incomplete,
                )

            argv = [
                "analyzer",
                "--local-eos", str(valid),
                "--additional-aa-local-eos", str(invalid),
                "--pp-local-eos", str(pp),
                "--out-dir", str(output),
                "--require-paired-aa",
                "--aa-events-per-chunk", "1",
            ]
            with patch.object(sys, "argv", argv), patch.object(
                analysis, "maybe_plot"
            ):
                self.assertEqual(analysis.main(), 0)

            with (output / "aa_sources.tsv").open() as handle:
                source_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [row["chunks_used"] for row in source_rows], ["1", "0"]
            )
            self.assertEqual(
                [row["chunks_skipped"] for row in source_rows], ["0", "1"]
            )
            with (output / "aa_rejections.tsv").open() as handle:
                rejection_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rejection_rows), 1)
            self.assertEqual(rejection_rows[0]["chunk_id"], "0")
            self.assertEqual(rejection_rows[0]["reason"], "ValueError")
            self.assertIn("no complete PYTHIA events", rejection_rows[0]["detail"])

    def test_v2_accepts_nonzero_hydro_slot_with_manifest_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            aa = root / "aa"
            pp = root / "pp"
            output = root / "analysis"
            self.make_snapshot(pp, "pp")
            (aa / "status/aa").mkdir(parents=True)
            (aa / "outputs/aa").mkdir(parents=True)
            (aa / "status/aa/chunk_0.txt").write_text("status=success\n")
            summary_header = (
                "variant\ttask_id\tseed\thydro_slot\thydro_event_id\t"
                "hydro_ncoll\thydro_payload_sha256\n"
            )
            digest = "a" * 64
            with tarfile.open(aa / "outputs/aa/chunk_0.tar.gz", "w:gz") as archive:
                for variant, suffix in (("no_prehydro", ""), ("with_prehydro", "_prehydro")):
                    parent = (
                        "runs/aa/hydro_026_C0-5_event_00643/"
                        f"task_00000{suffix}"
                    )
                    self.add_member(
                        archive, f"{parent}/HYBRID_Hadrons.out", self.event
                    )
                    summary = (
                        summary_header
                        + f"{variant}\t0\t900000\t26\t643\t35\t{digest}\n"
                    ).encode()
                    self.add_member(archive, f"{parent}/summary.tsv", summary)
            manifest = root / "aa_task_manifest.tsv"
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "task_id",
                        "hard_seed",
                        "hydro_slot",
                        "hydro_event_id",
                        "hydro_ncoll",
                        "hydro_payload_sha256",
                    ],
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "task_id": 0,
                        "hard_seed": 900000,
                        "hydro_slot": 26,
                        "hydro_event_id": 643,
                        "hydro_ncoll": 35,
                        "hydro_payload_sha256": digest,
                    }
                )
            argv = [
                "analyzer",
                "--local-eos", str(aa),
                "--pp-local-eos", str(pp),
                "--out-dir", str(output),
                "--require-paired-aa",
                "--aa-events-per-chunk", "1",
                "--aa-task-manifest", str(manifest),
                "--aa-task-limit", "1",
                "--require-complete-aa-prefix",
            ]
            with patch.object(sys, "argv", argv), patch.object(
                analysis, "maybe_plot"
            ):
                self.assertEqual(analysis.main(), 0)
            with (output / "aa_hydro_counts.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["hydro_slot"], "26")
            self.assertEqual(rows[0]["accepted_tasks"], "1")


if __name__ == "__main__":
    unittest.main()
