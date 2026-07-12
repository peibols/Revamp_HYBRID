#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import io
from pathlib import Path
import tarfile
import tempfile
import unittest


MODULE_PATH = Path(__file__).with_name("freeze_oo_v2_snapshot.py")
SPEC = importlib.util.spec_from_file_location("freeze_oo_v2_snapshot", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class FreezeSnapshotTest(unittest.TestCase):
    @staticmethod
    def add_text(archive: tarfile.TarFile, name: str, text: str) -> None:
        payload = text.encode()
        member = tarfile.TarInfo(name)
        member.size = len(payload)
        archive.addfile(member, io.BytesIO(payload))

    @classmethod
    def write_pair(cls, path: Path, task_id: int) -> None:
        task = f"task_{task_id:05d}"
        pair = (
            "kind\ttask_id\tvariant\treturncode\ttimeout\n"
            f"aa\t{task_id}\tno_prehydro\t0\t0\n"
            f"aa\t{task_id}\twith_prehydro\t0\t0\n"
        )
        summary = "variant\ttask_id\treturncode\ttimeout\n"
        with tarfile.open(path, "w:gz") as archive:
            cls.add_text(
                archive,
                f"run/{task}/summary.tsv",
                summary + f"no_prehydro\t{task_id}\t0\t0\n",
            )
            cls.add_text(archive, f"run/{task}/HYBRID_Hadrons.out", "event\n")
            cls.add_text(
                archive,
                f"run/{task}_prehydro/summary.tsv",
                summary + f"with_prehydro\t{task_id}\t0\t0\n",
            )
            cls.add_text(
                archive, f"run/{task}_prehydro/HYBRID_Hadrons.out", "event\n"
            )
            cls.add_text(archive, f"run/{task}_pair_summary.tsv", pair)

    def test_freezes_only_strict_ready_pairs_with_hardlinks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source"
            (source / "status/aa").mkdir(parents=True)
            (source / "outputs/aa").mkdir(parents=True)
            (source / "status/aa/chunk_1.txt").write_text("status=success\n")
            (source / "status/aa/chunk_2.txt").write_text("status=failed\n")
            self.write_pair(source / "outputs/aa/chunk_1.tar.gz", 1)
            self.write_pair(source / "outputs/aa/chunk_2.tar.gz", 2)
            manifest = root / "manifest.tsv"
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=("task_id", "hydro_slot"),
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(
                    ({"task_id": 1, "hydro_slot": 3}, {"task_id": 2, "hydro_slot": 4})
                )
            output = root / "frozen"
            metadata = MODULE.freeze_snapshot(
                source=source,
                manifest=manifest,
                output=output,
                source_repository=root,
            )
            self.assertEqual(metadata["strict_accepted_pairs"], 1)
            self.assertEqual((output / "accepted_task_ids.txt").read_text(), "1\n")
            selected = output / "local_eos/outputs/aa/chunk_1.tar.gz"
            self.assertTrue(selected.is_file())
            self.assertFalse(
                (output / "local_eos/outputs/aa/chunk_2.tar.gz").exists()
            )
            self.assertEqual(selected.stat().st_ino, (source / "outputs/aa/chunk_1.tar.gz").stat().st_ino)


if __name__ == "__main__":
    unittest.main()
