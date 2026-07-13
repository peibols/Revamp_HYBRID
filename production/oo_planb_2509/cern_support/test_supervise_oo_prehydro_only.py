#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import io
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest


SUPPORT = Path(__file__).parent
sys.path.insert(0, str(SUPPORT))
MODULE_PATH = SUPPORT / "supervise_oo_prehydro_only.py"
SPEC = importlib.util.spec_from_file_location("supervise_oo_prehydro_only", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class PrehydroOnlyArchiveTest(unittest.TestCase):
    assignment = {
        "task_id": "7",
        "hard_seed": "900007",
        "milestone_block": "1",
        "hydro_slot": "12",
        "hydro_event_id": "345",
        "hydro_ncoll": "2",
        "hydro_dir": "C0-5_event_00345",
        "hydro_payload_key": "hydro/C0-5/event_00345.tar.gz",
        "hydro_payload_sha256": "a" * 64,
    }

    @staticmethod
    def add_text(archive: tarfile.TarFile, name: str, text: str) -> None:
        payload = text.encode()
        member = tarfile.TarInfo(name)
        member.size = len(payload)
        archive.addfile(member, io.BytesIO(payload))

    def summary(self, alpha: float = 0.335) -> str:
        header = (
            "kind\ttask_id\tseed\tevents\tcentrality\thydro_index\t"
            "hydro_event_id\thydro_ncoll\thydro_payload_sha256\tvariant\t"
            "use_prehydro\tenergy_loss_alpha\tbroadening_k\tprehydro_file\t"
            "returncode\ttimeout\tseconds\tdir\n"
        )
        row = (
            f"aa\t7\t900007\t1\tC0-5\t12\t345\t2\t{'a' * 64}\t"
            f"with_prehydro\t1\t{alpha}\t15.0\tprehydro_table.tsv\t"
            "0\t0\t1.0\truns/test/task_00007_prehydro\n"
        )
        return header + row

    def make_archive(
        self,
        path: Path,
        *,
        alpha: float = 0.335,
        include_baseline: bool = False,
    ) -> None:
        base = "runs/test/aa/hydro"
        with tarfile.open(path, "w:gz") as archive:
            self.add_text(
                archive,
                f"{base}/task_00007_prehydro/summary.tsv",
                self.summary(alpha),
            )
            self.add_text(
                archive,
                f"{base}/task_00007_prehydro/HYBRID_Hadrons.out",
                "# event 0\nweight 1 cross 1\nend\n",
            )
            self.add_text(
                archive,
                f"{base}/task_00007_prehydro/hybrid_input.dat",
                f"alpha = {alpha}\nuse_prehydro = true\n",
            )
            self.add_text(
                archive,
                f"{base}/task_00007_prehydro_only_summary.tsv",
                self.summary(alpha),
            )
            if include_baseline:
                self.add_text(
                    archive,
                    f"{base}/task_00007/HYBRID_Hadrons.out",
                    "unexpected baseline\n",
                )

    def test_accepts_exact_single_variant_archive(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chunk_7.tar.gz"
            self.make_archive(path)
            self.assertTrue(
                MODULE.prehydro_only_archive_is_complete(
                    path,
                    self.assignment,
                    expected_alpha=0.335,
                    expected_broadening_k=15.0,
                )
            )

    def test_rejects_wrong_alpha(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chunk_7.tar.gz"
            self.make_archive(path, alpha=0.37)
            self.assertFalse(
                MODULE.prehydro_only_archive_is_complete(
                    path,
                    self.assignment,
                    expected_alpha=0.335,
                    expected_broadening_k=15.0,
                )
            )

    def test_rejects_archive_containing_baseline_leg(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chunk_7.tar.gz"
            self.make_archive(path, include_baseline=True)
            self.assertFalse(
                MODULE.prehydro_only_archive_is_complete(
                    path,
                    self.assignment,
                    expected_alpha=0.335,
                    expected_broadening_k=15.0,
                )
            )

    def test_status_requires_prehydro_only_provenance(self) -> None:
        status = {
            "kind": "aa",
            "task_id": "7",
            "status": "success",
            "exit_code": "0",
            "run_prehydro_pair": "false",
            "run_prehydro_only": "true",
            "prehydro_alpha": "0.335",
            "broadening_k": "15.0",
            "hydro_slot": "12",
            "hydro_event_id": "345",
            "hydro_ncoll": "2",
            "hydro_dir": "C0-5_event_00345",
            "hydro_payload_key": "hydro/C0-5/event_00345.tar.gz",
            "hydro_payload_sha256": "a" * 64,
        }
        self.assertTrue(
            MODULE.status_matches_assignment(
                status,
                self.assignment,
                expected_alpha=0.335,
                expected_broadening_k=15.0,
            )
        )
        status["run_prehydro_only"] = "false"
        self.assertFalse(
            MODULE.status_matches_assignment(
                status,
                self.assignment,
                expected_alpha=0.335,
                expected_broadening_k=15.0,
            )
        )


if __name__ == "__main__":
    unittest.main()
