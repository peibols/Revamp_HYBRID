#!/usr/bin/env python3

from __future__ import annotations

import io
from pathlib import Path
import tarfile
import tempfile
import unittest

import supervise_oo_moliere_pair as supervisor


TABLE_SHA = "3" * 64


class MolierePairAuditTest(unittest.TestCase):
    assignment = {
        "task_id": "7",
        "hard_seed": "900007",
        "hydro_slot": "12",
        "hydro_event_id": "345",
        "hydro_ncoll": "3",
        "hydro_dir": "C0-5_event_00345",
        "hydro_payload_key": "hydro/C0-5/event_00345.tar.gz",
        "hydro_payload_sha256": "a" * 64,
    }

    @staticmethod
    def add_bytes(archive: tarfile.TarFile, name: str, payload: bytes) -> None:
        member = tarfile.TarInfo(name)
        member.size = len(payload)
        archive.addfile(member, io.BytesIO(payload))

    def row(self, variant: str, use_prehydro: bool) -> bytes:
        header = (
            "variant\tuse_prehydro\tenergy_loss_alpha\tbroadening_k\t"
            "do_moliere\thadro_type\tmoliere_tables_sha256\tprehydro_file\t"
            "returncode\ttimeout\tseconds\tdir\ttask_id\tseed\tcentrality\t"
            "hydro_slot\thydro_event_id\thydro_ncoll\thydro_payload_sha256\n"
        )
        payload = (
            f"{variant}\t{int(use_prehydro)}\t0.335\t15.0\t1\t1\t{TABLE_SHA}\t"
            f"{'prehydro_table.tsv' if use_prehydro else ''}\t0\t0\t1.0\trun\t"
            f"7\t900007\tC0-5\t12\t345\t3\t{'a' * 64}\n"
        )
        return (header + payload).encode()

    @staticmethod
    def input_card(*, use_prehydro: bool, do_elastic: bool = True) -> bytes:
        values = {
            "alpha": "0.335",
            "do_quench": "true",
            "do_wake": "true",
            "do_elastic": "true" if do_elastic else "false",
            "do_lres": "false",
            "do_Moliere_on_unresolved_partons": "false",
            "do_Moliere_dynamic_unresolved_resolution": "false",
            "do_Moliere_dynamic_daughter_unresolved_resolution": "false",
            "do_Moliere_recursive_unresolved_resolution": "false",
            "hadro_type": "1",
            "compat_moliere_legacy_hydro": "true",
            "use_prehydro": "true" if use_prehydro else "false",
            "tables_path": "/work/a10_tables/",
        }
        return "".join(f"{key} = {value}\n" for key, value in values.items()).encode()

    def make_archive(self, path: Path, *, do_elastic: bool = True) -> None:
        base = "runs/test/aa/hydro/task_00007"
        pre = f"{base}_prehydro"
        with tarfile.open(path, "w:gz") as archive:
            self.add_bytes(archive, f"{base}/summary.tsv", self.row("no_prehydro", False))
            self.add_bytes(archive, f"{base}/HYBRID_Hadrons.out", b"hadron\n")
            self.add_bytes(
                archive,
                f"{base}/hybrid_input.dat",
                self.input_card(use_prehydro=False, do_elastic=do_elastic),
            )
            self.add_bytes(archive, f"{pre}/summary.tsv", self.row("with_prehydro", True))
            self.add_bytes(archive, f"{pre}/HYBRID_Hadrons.out", b"hadron\n")
            self.add_bytes(
                archive,
                f"{pre}/hybrid_input.dat",
                self.input_card(use_prehydro=True, do_elastic=do_elastic),
            )
            pair = self.row("no_prehydro", False) + self.row("with_prehydro", True).split(b"\n", 1)[1]
            self.add_bytes(archive, f"runs/test/aa/hydro/task_00007_pair_summary.tsv", pair)

    def test_accepts_strict_moliere_pair(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chunk_7.tar.gz"
            self.make_archive(path)
            self.assertTrue(
                supervisor.archive_is_complete(
                    path,
                    self.assignment,
                    alpha=0.335,
                    broadening_k=15.0,
                    tables_sha256=TABLE_SHA,
                )
            )

    def test_rejects_archive_with_elastic_disabled(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "chunk_7.tar.gz"
            self.make_archive(path, do_elastic=False)
            self.assertFalse(
                supervisor.archive_is_complete(
                    path,
                    self.assignment,
                    alpha=0.335,
                    broadening_k=15.0,
                    tables_sha256=TABLE_SHA,
                )
            )

    def test_status_requires_moliere_provenance(self) -> None:
        status = {
            "kind": "aa",
            "task_id": "7",
            "status": "success",
            "exit_code": "0",
            "run_prehydro_pair": "true",
            "run_prehydro_only": "false",
            "no_prehydro_alpha": "0.335",
            "prehydro_alpha": "0.335",
            "broadening_k": "15.0",
            "do_moliere": "true",
            "moliere_mode": "legacy_resolved",
            "moliere_tables_sha256": TABLE_SHA,
            "hydro_slot": "12",
            "hydro_event_id": "345",
            "hydro_ncoll": "3",
            "hydro_dir": "C0-5_event_00345",
            "hydro_payload_key": "hydro/C0-5/event_00345.tar.gz",
            "hydro_payload_sha256": "a" * 64,
        }
        self.assertTrue(
            supervisor.status_matches(
                status,
                self.assignment,
                alpha=0.335,
                broadening_k=15.0,
                tables_sha256=TABLE_SHA,
            )
        )
        status["do_moliere"] = "false"
        self.assertFalse(
            supervisor.status_matches(
                status,
                self.assignment,
                alpha=0.335,
                broadening_k=15.0,
                tables_sha256=TABLE_SHA,
            )
        )


if __name__ == "__main__":
    unittest.main()
