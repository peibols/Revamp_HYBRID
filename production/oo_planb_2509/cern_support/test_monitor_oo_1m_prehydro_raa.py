#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import csv
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

    def test_accepts_live_additional_source_sync(self) -> None:
        argv = [
            "monitor",
            "--additional-aa-local-eos", "/tmp/aa20k",
            "--additional-aa-chunks", "20000",
            "--sync-via-cernctl",
            "--sync-additional-aa", "/eos/aa20k", "/tmp/aa20k",
        ]
        with patch.object(sys, "argv", argv):
            args = monitor.parse_args()
        self.assertEqual(
            args.sync_additional_aa,
            [("/eos/aa20k", Path("/tmp/aa20k"))],
        )

    def test_rejects_sync_path_not_in_additional_sources(self) -> None:
        argv = [
            "monitor",
            "--sync-via-cernctl",
            "--sync-additional-aa", "/eos/aa20k", "/tmp/aa20k",
        ]
        with patch.object(sys, "argv", argv), self.assertRaises(SystemExit):
            monitor.parse_args()

    def test_syncs_primary_and_live_additional_sources(self) -> None:
        args = SimpleNamespace(
            cernctl="cernctl",
            cern_remote="lxplus",
            remote_stage="/tmp/snapshot.tar.gz",
            renew_kerberos=True,
            sync_additional_aa=[("/eos/aa20k", Path("/tmp/aa20k"))],
        )
        with patch.object(monitor, "sync_via_cernctl") as sync:
            monitor.sync_configured_aa_sources(
                args,
                eos_base="/eos/primary",
                local_eos=Path("/tmp/primary"),
                include_outputs=True,
            )
        self.assertEqual(sync.call_count, 2)
        primary_call, additional_call = sync.call_args_list
        self.assertTrue(primary_call.kwargs["renew_kerberos"])
        self.assertFalse(additional_call.kwargs["renew_kerberos"])
        self.assertEqual(additional_call.kwargs["eos_base"], "/eos/aa20k")
        self.assertEqual(additional_call.kwargs["local_eos"], Path("/tmp/aa20k"))

    def test_live_sync_uses_incremental_rsync(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            local_eos = Path(tmp) / "snapshot"
            with patch.object(monitor.subprocess, "run") as run:
                monitor.sync_via_cernctl(
                    cernctl="cernctl",
                    cern_remote="lxplus",
                    eos_base="/eos/aa20k",
                    local_eos=local_eos,
                    remote_stage="/tmp/legacy.tar.gz",
                    include_outputs=True,
                    renew_kerberos=True,
                )
        commands = [call.args[0] for call in run.call_args_list]
        self.assertEqual(commands[0], ["cernctl", "kerberos"])
        self.assertEqual(
            commands[1],
            [
                "rsync", "-a", "--partial",
                "lxplus:/eos/aa20k/status/",
                f"{local_eos / 'status'}/",
            ],
        )
        self.assertEqual(
            commands[2],
            [
                "rsync", "-a", "--partial",
                "lxplus:/eos/aa20k/outputs/",
                f"{local_eos / 'outputs'}/",
            ],
        )

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

    def test_reads_only_analyzer_accepted_aa_chunks(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            out_dir = Path(tmp)
            with (out_dir / "aa_sources.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(
                    ["aa_local_eos", "chunks_used", "chunks_skipped", "missing_outputs"]
                )
                writer.writerow(["/tmp/aa10k", 9999, 1, 0])
                writer.writerow(["/tmp/aa20k", 16919, 1, 0])
            self.assertEqual(monitor.read_analyzed_aa_chunks(out_dir), 26918)

    def test_selects_only_highest_newly_reached_milestone(self) -> None:
        milestones = [10, 20, 30, 40, 50, 60]
        self.assertEqual(
            monitor.highest_pending_milestone(milestones, 56, set()),
            50,
        )
        self.assertEqual(
            monitor.highest_pending_milestone(milestones, 56, {50}),
            None,
        )
        self.assertEqual(
            monitor.highest_pending_milestone(milestones, 64, {50}),
            60,
        )


if __name__ == "__main__":
    unittest.main()
