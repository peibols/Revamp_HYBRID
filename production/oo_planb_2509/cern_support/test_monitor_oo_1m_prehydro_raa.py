#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import monitor_oo_1m_prehydro_raa as monitor


class CombinedMonitorArgumentsTest(unittest.TestCase):
    def test_accepts_matching_additional_sources_and_counts(self) -> None:
        argv = [
            "monitor",
            "--additional-aa-local-eos", "/tmp/aa10k",
            "--additional-aa-chunks", "10000",
        ]
        with patch.object(sys, "argv", argv):
            args = monitor.parse_args()
        self.assertEqual(args.additional_aa_local_eos, [Path("/tmp/aa10k")])
        self.assertEqual(args.additional_aa_chunks, [10000])

    def test_rejects_unmatched_additional_source(self) -> None:
        argv = ["monitor", "--additional-aa-local-eos", "/tmp/aa10k"]
        with patch.object(sys, "argv", argv), self.assertRaises(SystemExit):
            monitor.parse_args()

    def test_forwards_all_aa_sources_to_analyzer(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            with patch.object(monitor.subprocess, "run") as run:
                monitor.run_analyzer(
                    work,
                    Path("/tmp/aa20k"),
                    [Path("/tmp/aa10k")],
                    Path("/tmp/pp1m"),
                    Path("/tmp/out"),
                    1,
                )
            command = run.call_args.args[0]
        self.assertIn("--additional-aa-local-eos", command)
        index = command.index("--additional-aa-local-eos")
        self.assertEqual(command[index + 1], "/tmp/aa10k")


if __name__ == "__main__":
    unittest.main()
