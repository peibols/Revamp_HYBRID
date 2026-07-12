#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import io
from pathlib import Path
import tarfile
import tempfile
import unittest
from unittest.mock import patch


MODULE_PATH = Path(__file__).with_name("supervise_oo_v2.py")
SPEC = importlib.util.spec_from_file_location("supervise_oo_v2", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class V2SupervisorTest(unittest.TestCase):
    @staticmethod
    def write_paired_archive(path: Path, task_id: int) -> None:
        task = f"task_{task_id:05d}"
        pair_header = (
            "kind\ttask_id\tvariant\treturncode\ttimeout\n"
            f"aa\t{task_id}\tno_prehydro\t0\t0\n"
            f"aa\t{task_id}\twith_prehydro\t0\t0\n"
        )
        variant_header = "variant\ttask_id\treturncode\ttimeout\n"
        files = {
            f"run/{task}/summary.tsv": (
                variant_header + f"no_prehydro\t{task_id}\t0\t0\n"
            ),
            f"run/{task}/HYBRID_Hadrons.out": "event\n",
            f"run/{task}_prehydro/summary.tsv": (
                variant_header + f"with_prehydro\t{task_id}\t0\t0\n"
            ),
            f"run/{task}_prehydro/HYBRID_Hadrons.out": "event\n",
            f"run/{task}_pair_summary.tsv": pair_header,
        }
        with tarfile.open(path, "w:gz") as archive:
            for name, text in files.items():
                payload = text.encode()
                member = tarfile.TarInfo(name)
                member.size = len(payload)
                archive.addfile(member, io.BytesIO(payload))

    def test_ready_ids_require_success_status_and_nonempty_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "status/aa").mkdir(parents=True)
            (root / "outputs/aa").mkdir(parents=True)
            (root / "status/aa/chunk_1.txt").write_text("status=success\n")
            (root / "status/aa/chunk_2.txt").write_text("status=failed\n")
            self.write_paired_archive(root / "outputs/aa/chunk_1.tar.gz", 1)
            self.write_paired_archive(root / "outputs/aa/chunk_2.tar.gz", 2)
            self.write_paired_archive(root / "outputs/aa/chunk_3.tar.gz", 3)
            self.assertEqual(MODULE.ready_ids(root), {1})

    def test_partial_pair_archive_is_not_ready(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "status/aa").mkdir(parents=True)
            (root / "outputs/aa").mkdir(parents=True)
            (root / "status/aa/chunk_4.txt").write_text("status=success\n")
            with tarfile.open(root / "outputs/aa/chunk_4.tar.gz", "w:gz") as archive:
                payload = b"event\n"
                member = tarfile.TarInfo("run/task_00004/HYBRID_Hadrons.out")
                member.size = len(payload)
                archive.addfile(member, io.BytesIO(payload))
            self.assertEqual(MODULE.ready_ids(root), set())

    def test_retry_submit_replaces_queue_and_log_names(self) -> None:
        template = """output = log/aa.$(ClusterId).out
error = log/aa.$(ClusterId).err
log = log/aa.$(ClusterId).log
when_to_transfer_output = ON_EXIT
queue chunk_id from aa_chunk_ids.txt
"""
        rendered = MODULE.render_retry_submit(
            template + 'environment = "TIMEOUT_S=39600"\n',
            "retry_ids.txt",
            "stamp",
            retry_timeout_s=72_000,
        )
        self.assertIn("queue chunk_id from retry_ids.txt", rendered)
        self.assertIn("output = log/v2_retry_stamp_aa.$(ClusterId).out", rendered)
        self.assertIn("error = log/v2_retry_stamp_aa.$(ClusterId).err", rendered)
        self.assertIn("log = /dev/null", rendered)
        self.assertIn("TIMEOUT_S=72000", rendered)
        self.assertIn('transfer_output_files = ""', rendered)

    def test_batched_ids_are_sorted_and_respect_submission_limit(self) -> None:
        batches = MODULE.batched_ids({8, 1, 5, 2, 9}, 2)
        self.assertEqual(batches, [[1, 2], [5, 8], [9]])

    def test_retry_submit_accepts_discarded_stdio(self) -> None:
        template = """output = /dev/null
error = /dev/null
log = log/aa.$(ClusterId).log
when_to_transfer_output = ON_EXIT
queue chunk_id from aa_chunk_ids.txt
"""
        rendered = MODULE.render_retry_submit(
            template + 'environment = "TIMEOUT_S=39600"\n',
            "retry_ids.txt",
            "stamp",
            retry_timeout_s=72_000,
        )
        self.assertIn("output = /dev/null", rendered)
        self.assertIn("error = /dev/null", rendered)
        self.assertIn("log = /dev/null", rendered)
        self.assertIn("queue chunk_id from retry_ids.txt", rendered)
        self.assertIn('transfer_output_files = ""', rendered)

    def test_dry_run_reports_but_retains_terminal_holds(self) -> None:
        args = type("Args", (), {"dry_run": True})()
        self.assertEqual(
            MODULE.remove_terminal_holds(args, {3, 7}, {3}),
            set(),
        )

    def test_terminal_holds_are_removed_for_retry_or_cleanup(self) -> None:
        args = type(
            "Args",
            (),
            {
                "dry_run": False,
                "campaign": "campaign",
                "schedd": "schedd",
                "cernctl": "cernctl",
            },
        )()
        with patch.object(MODULE.subprocess, "run") as run:
            removed = MODULE.remove_terminal_holds(args, {3, 7}, {3})
        self.assertEqual(removed, {3, 7})
        run.assert_called_once()
        command = run.call_args.args[0]
        self.assertEqual(command[:3], ["cernctl", "run", "bash"])
        self.assertIn("condor_rm", command[-1])

    def test_shifted_analysis_passes_global_task_range(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            args = type(
                "Args",
                (),
                {
                    "work": root / "work",
                    "source": root / "source",
                    "pp_local_eos": root / "pp",
                    "task_id_start": 50_000,
                },
            )()
            with patch.object(MODULE.subprocess, "run") as run:
                output = MODULE.run_analysis(args, root / "eos", 55_000)
            command = run.call_args.args[0]
            self.assertEqual(
                output,
                root / "work/analysis_v2_milestones/tasks_50000_55000",
            )
            self.assertEqual(
                command[command.index("--aa-task-start") + 1], "50000"
            )
            self.assertEqual(
                command[command.index("--aa-task-limit") + 1], "55000"
            )


if __name__ == "__main__":
    unittest.main()
