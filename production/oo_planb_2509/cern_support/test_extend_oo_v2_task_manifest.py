#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


MODULE_PATH = Path(__file__).with_name("extend_oo_v2_task_manifest.py")
SPEC = importlib.util.spec_from_file_location("extend_oo_v2_task_manifest", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class ExtendManifestTest(unittest.TestCase):
    def test_shift_preserves_hydro_and_separates_task_and_seed_ranges(self) -> None:
        rows = [
            {
                "task_id": str(task_id),
                "hard_seed": str(900_000 + task_id),
                "milestone_block": "1",
                "hydro_slot": str(task_id % 2),
                "hydro_event_id": str(10 + task_id % 2),
                "hydro_ncoll": "20",
            }
            for task_id in range(4)
        ]
        shifted = MODULE.shifted_rows(
            rows,
            task_id_offset=50_000,
            milestone_block_offset=10,
        )
        self.assertEqual(
            [int(row["task_id"]) for row in shifted],
            list(range(50_000, 50_004)),
        )
        self.assertEqual(
            [int(row["hard_seed"]) for row in shifted],
            list(range(950_000, 950_004)),
        )
        self.assertTrue(all(row["milestone_block"] == "11" for row in shifted))
        self.assertEqual(
            [row["hydro_event_id"] for row in shifted],
            [row["hydro_event_id"] for row in rows],
        )

    def test_reader_rejects_noncontiguous_task_ids(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.tsv"
            fieldnames = [*MODULE.INTEGER_FIELDS]
            with path.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=fieldnames,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerow(dict(zip(fieldnames, (0, 10, 1, 0, 1, 2))))
                writer.writerow(dict(zip(fieldnames, (2, 11, 1, 0, 1, 2))))
            with self.assertRaisesRegex(ValueError, "ordered and contiguous"):
                MODULE.read_manifest(path)


if __name__ == "__main__":
    unittest.main()
