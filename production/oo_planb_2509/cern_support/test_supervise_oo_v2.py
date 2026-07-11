#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest


MODULE_PATH = Path(__file__).with_name("supervise_oo_v2.py")
SPEC = importlib.util.spec_from_file_location("supervise_oo_v2", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class V2SupervisorTest(unittest.TestCase):
    def test_ready_ids_require_success_status_and_nonempty_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "status/aa").mkdir(parents=True)
            (root / "outputs/aa").mkdir(parents=True)
            (root / "status/aa/chunk_1.txt").write_text("status=success\n")
            (root / "status/aa/chunk_2.txt").write_text("status=failed\n")
            (root / "outputs/aa/chunk_1.tar.gz").write_bytes(b"output")
            (root / "outputs/aa/chunk_2.tar.gz").write_bytes(b"output")
            (root / "outputs/aa/chunk_3.tar.gz").write_bytes(b"output")
            self.assertEqual(MODULE.ready_ids(root), {1})

    def test_retry_submit_replaces_queue_and_log_names(self) -> None:
        template = """output = log/aa.$(ClusterId).out
error = log/aa.$(ClusterId).err
log = log/aa.$(ClusterId).log
queue chunk_id from aa_chunk_ids.txt
"""
        rendered = MODULE.render_retry_submit(template, "retry_ids.txt", "stamp")
        self.assertIn("queue chunk_id from retry_ids.txt", rendered)
        self.assertIn("output = log/v2_retry_stamp_aa.$(ClusterId).out", rendered)
        self.assertIn("error = log/v2_retry_stamp_aa.$(ClusterId).err", rendered)


if __name__ == "__main__":
    unittest.main()
