#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import supervise_oo_50k as supervisor


class ReconciliationHelpersTest(unittest.TestCase):
    def test_reads_union_of_expected_id_files(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = root / "first.txt"
            second = root / "second.txt"
            first.write_text("1\n2\n")
            second.write_text("2\n3\n")
            self.assertEqual(supervisor.read_ids((first, second)), {1, 2, 3})

    def test_ready_ids_require_success_and_nonempty_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            local_eos = Path(tmp)
            status = local_eos / "status/aa"
            outputs = local_eos / "outputs/aa"
            status.mkdir(parents=True)
            outputs.mkdir(parents=True)
            (status / "chunk_1.txt").write_text("status=success\n")
            (status / "chunk_2.txt").write_text("status=failed\n")
            (status / "chunk_3.txt").write_text("status=success\n")
            (outputs / "chunk_1.tar.gz").write_bytes(b"archive")
            (outputs / "chunk_2.tar.gz").write_bytes(b"archive")
            (outputs / "chunk_3.tar.gz").write_bytes(b"")
            self.assertEqual(supervisor.ready_ids(local_eos), {1})

    def test_retry_submit_uses_new_id_file_and_unique_logs(self) -> None:
        template = """output = log/run.$(ClusterId).out
error = log/run.$(ClusterId).err
log = log/run.$(ClusterId).log
queue chunk_id from old_ids.txt
"""
        rendered = supervisor.render_retry_submit(
            template, "retry_ids.txt", "20260711T100000"
        )
        self.assertIn("queue chunk_id from retry_ids.txt", rendered)
        self.assertNotIn("old_ids.txt", rendered)
        self.assertEqual(rendered.count("reconcile_20260711T100000_"), 3)


if __name__ == "__main__":
    unittest.main()
