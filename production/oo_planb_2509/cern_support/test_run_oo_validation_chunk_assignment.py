#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest


MODULE_PATH = Path(__file__).with_name("run_oo_validation_chunk.py")
SPEC = importlib.util.spec_from_file_location("run_oo_validation_chunk", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class V2AssignmentTest(unittest.TestCase):
    def make_manifest(self, path: Path, **updates: object) -> None:
        row = {
            "task_id": 7,
            "hard_seed": 900007,
            "milestone_block": 1,
            "hydro_slot": 12,
            "hydro_event_id": 345,
            "hydro_ncoll": 3,
            "hydro_dir": "C0-5_event_00345",
            "hydro_payload_key": "hydro/C0-5/event_00345.tar.gz",
            "hydro_payload_sha256": "a" * 64,
        }
        row.update(updates)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=list(row), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerow(row)

    def test_manifest_and_staged_hydro_match(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "tasks.tsv"
            self.make_manifest(manifest)
            assignment = MODULE.load_aa_task_assignment(manifest, 7)
            self.assertEqual(assignment["hard_seed"], 900007)
            self.assertEqual(assignment["hydro_event_id"], 345)

            hydro = root / "C0-5_event_00345"
            hydro.mkdir()
            (hydro / "README_staged_event.txt").write_text(
                "hydro_slot = 12\nevent_id = 345\nncoll_positions = 3\n"
            )
            (hydro / "NcollList.dat").write_text("# x y\n0 0\n1 1\n2 2\n")
            self.assertEqual(
                MODULE.read_staged_hydro_metadata(hydro),
                {"hydro_slot": 12, "event_id": 345, "ncoll_positions": 3},
            )

    def test_manifest_rejects_path_traversal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = Path(tmp) / "tasks.tsv"
            self.make_manifest(manifest, hydro_dir="../C0-5_event_00345")
            with self.assertRaisesRegex(ValueError, "hydro_dir"):
                MODULE.load_aa_task_assignment(manifest, 7)

    def test_staged_hydro_rejects_ncoll_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            hydro = Path(tmp)
            (hydro / "README_staged_event.txt").write_text(
                "hydro_slot = 12\nevent_id = 345\nncoll_positions = 3\n"
            )
            (hydro / "NcollList.dat").write_text("0 0\n1 1\n")
            with self.assertRaisesRegex(ValueError, "contains 2 positions"):
                MODULE.read_staged_hydro_metadata(hydro)

    def test_write_input_records_requested_energy_loss_parameters(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "hybrid_input.dat"
            MODULE.write_input(
                path,
                seed=900007,
                events=1,
                aa=True,
                energy_loss_alpha=0.355,
                broadening_k=15.0,
                use_prehydro=True,
            )
            values = {
                key.strip(): value.strip()
                for key, value in (
                    line.split("=", 1)
                    for line in path.read_text().splitlines()
                    if "=" in line
                )
            }
            self.assertEqual(values["alpha"], "0.355")
            self.assertEqual(values["kappa"], "15.0")
            self.assertEqual(values["use_prehydro"], "true")


if __name__ == "__main__":
    unittest.main()
