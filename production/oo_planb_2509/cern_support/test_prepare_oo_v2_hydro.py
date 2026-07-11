#!/usr/bin/env python3

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import io
from pathlib import Path
import struct
import tarfile
import tempfile
import unittest
from unittest.mock import patch


MODULE_PATH = Path(__file__).with_name("prepare_oo_v2_hydro.py")
SPEC = importlib.util.spec_from_file_location("prepare_oo_v2_hydro", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def hydro_payload() -> bytes:
    values = [0.4, 0.1, 150.0, 0.2, -15.0, 150.0, 0.2, -15.0]
    values.extend([0.0] * 7)
    values.append(11.0)
    return struct.pack("<16f", *values) + b"payload"


def add_bytes(archive: tarfile.TarFile, name: str, payload: bytes) -> None:
    member = tarfile.TarInfo(name)
    member.size = len(payload)
    archive.addfile(member, io.BytesIO(payload))


class HydroPreparationTest(unittest.TestCase):
    def make_archive(self, path: Path) -> str:
        with tarfile.open(path, "w:gz") as archive:
            for event_id, ncoll in ((10, 2), (20, 3)):
                base = f"C0-5/hydro_results_{event_id}"
                add_bytes(archive, f"{base}/evolution_all_xyeta.dat", hydro_payload())
                positions = b"# x y\n" + b"0 0\n" * ncoll
                add_bytes(archive, f"{base}/NcollList{event_id}.dat", positions)
                add_bytes(archive, f"{base}/music_input", b"music\n")
                add_bytes(archive, f"{base}/run.log", b"run\n")
        return hashlib.md5(path.read_bytes()).hexdigest()

    def test_largest_remainder_closes(self) -> None:
        self.assertEqual(MODULE.largest_remainder_quotas([2, 3], 10), [4, 6])

    def test_end_to_end_payload_and_balanced_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "C0-5.tar.gz"
            expected_md5 = self.make_archive(source)
            output = root / "output"
            argv = [
                "prepare",
                "--archive", str(source),
                "--out-dir", str(output),
                "--expected-events", "2",
                "--expected-md5", expected_md5,
                "--hard-events", "20",
                "--milestone-block-size", "10",
                "--seed-offset", "100",
                "--shuffle-seed", "7",
            ]
            with patch("sys.argv", argv):
                self.assertEqual(MODULE.main(), 0)

            with (output / "hydro_manifest.tsv").open() as handle:
                hydros = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["event_id"] for row in hydros], ["10", "20"])
            self.assertEqual([row["ncoll"] for row in hydros], ["2", "3"])

            with (output / "aa_task_manifest.tsv").open() as handle:
                tasks = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(tasks), 20)
            self.assertEqual([int(row["hard_seed"]) for row in tasks], list(range(100, 120)))
            for block in (1, 2):
                selected = [row for row in tasks if int(row["milestone_block"]) == block]
                counts = {event_id: 0 for event_id in (10, 20)}
                for row in selected:
                    counts[int(row["hydro_event_id"])] += 1
                self.assertEqual(counts, {10: 4, 20: 6})

            payload = output / "hydro_payloads/event_00010.tar.gz"
            with gzip.open(payload, "rb") as stream:
                with tarfile.open(fileobj=stream, mode="r:") as archive:
                    names = set(archive.getnames())
            self.assertIn("C0-5_event_00010/evolution_all_xyeta.dat", names)
            self.assertIn("C0-5_event_00010/NcollList.dat", names)


if __name__ == "__main__":
    unittest.main()
