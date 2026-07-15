#!/usr/bin/env python3

from __future__ import annotations

import hashlib
import io
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest
import zipfile


RUNNER = Path(__file__).with_name("run_chunk_job.sh")


class ChunkFailureToleranceTest(unittest.TestCase):
    def run_fetch_failure(self, tolerate: bool) -> tuple[subprocess.CompletedProcess[str], str]:
        with tempfile.TemporaryDirectory() as tmp:
            work = Path(tmp)
            fake_bin = work / "bin"
            fake_bin.mkdir()
            xrdcp = fake_bin / "xrdcp"
            xrdcp.write_text("#!/usr/bin/env bash\nexit 1\n")
            xrdcp.chmod(0o755)
            env = os.environ.copy()
            env.update({
                "PATH": f"{fake_bin}:{env['PATH']}",
                "EOS_BASE": "/eos/test/campaign",
                "KIND": "aa",
                "SEED_OFFSET": "830000",
                "EVENTS": "1",
                "RUN_NAME": "runs/test",
                "RUN_PREHYDRO_PAIR": "false",
                "TOLERATE_CHUNK_FAILURE": str(tolerate).lower(),
            })
            result = subprocess.run(
                ["bash", str(RUNNER), "7"],
                cwd=work,
                env=env,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            status = (work / "chunk_status.txt").read_text()
        return result, status

    def test_default_behavior_propagates_chunk_failure(self) -> None:
        result, status = self.run_fetch_failure(tolerate=False)
        self.assertEqual(result.returncode, 1)
        self.assertIn("status=failed", status)
        self.assertIn("exit_code=1", status)

    def test_tolerated_failure_is_retained_but_does_not_fail_condor(self) -> None:
        result, status = self.run_fetch_failure(tolerate=True)
        self.assertEqual(result.returncode, 0)
        self.assertIn("status=failed", status)
        self.assertIn("exit_code=1", status)
        self.assertIn("failure retained in EOS status", result.stderr)

    @staticmethod
    def add_bytes(archive: tarfile.TarFile, name: str, payload: bytes) -> None:
        member = tarfile.TarInfo(name)
        member.size = len(payload)
        archive.addfile(member, io.BytesIO(payload))

    def run_v2_wrapper(
        self,
        valid_checksum: bool,
        *,
        prehydro_only: bool = False,
        do_moliere: bool = False,
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            eos = root / "mock_eos/eos/test/campaign"
            payloads = root / "mock_eos/eos/test/shared/payloads"
            runtime_payloads = root / "mock_eos/eos/test/runtime/payloads"
            moliere_payloads = root / "mock_eos/eos/test/moliere/payloads"
            (payloads / "hydro/C0-5").mkdir(parents=True)
            runtime_payloads.mkdir(parents=True)
            moliere_payloads.mkdir(parents=True)
            fake_bin = root / "fake_bin"
            fake_bin.mkdir()
            xrdcp = fake_bin / "xrdcp"
            xrdcp.write_text(
                """#!/usr/bin/env bash
set -euo pipefail
src="${@: -2:1}"
dst="${@: -1}"
map_path() {
  local value="$1"
  if [[ "$value" == root://eosuser.cern.ch/* ]]; then
    value="${value#root://eosuser.cern.ch/}"
    value="${value#/}"
    printf '%s/%s' "$MOCK_EOS" "$value"
  else
    printf '%s' "$value"
  fi
}
source_path="$(map_path "$src")"
destination_path="$(map_path "$dst")"
mkdir -p "$(dirname "$destination_path")"
cp "$source_path" "$destination_path"
"""
            )
            xrdcp.chmod(0o755)

            pythia_root = root / "pythia8315_alma9_install"
            (pythia_root / "share/Pythia8/xmldoc").mkdir(parents=True)
            (pythia_root / "lib").mkdir()
            with tarfile.open(payloads / "pythia8315_alma9_install.tar.gz", "w:gz") as archive:
                archive.add(pythia_root, arcname=pythia_root.name)

            runtime_root = root / "runtime_build"
            (runtime_root / "runtime").mkdir(parents=True)
            (runtime_root / "bin").mkdir()
            fake_runner = runtime_root / "runtime/run_oo_validation_chunk.py"
            fake_runner.write_text(
                """#!/usr/bin/env python3
import argparse
from pathlib import Path
p = argparse.ArgumentParser(add_help=False)
p.add_argument('--run-name')
p.add_argument('--aa-task-manifest')
p.add_argument('--no-prehydro-alpha')
p.add_argument('--prehydro-alpha')
p.add_argument('--broadening-k')
p.add_argument('--run-prehydro-only', action='store_true')
p.add_argument('--do-moliere', action='store_true')
p.add_argument('--moliere-tables-path')
p.add_argument('--moliere-tables-sha256')
args, _ = p.parse_known_args()
assert Path(args.aa_task_manifest).is_file()
assert args.no_prehydro_alpha == '0.37'
assert args.prehydro_alpha == '0.355'
assert args.broadening_k == '15.0'
assert args.run_prehydro_only == PREHYDRO_ONLY
assert args.do_moliere == DO_MOLIERE
if args.do_moliere:
    table_root = Path(args.moliere_tables_path)
    assert len(list((table_root / 'quark_tables').glob('*.dat'))) == 476
    assert len(list((table_root / 'gluon_tables').glob('*.dat'))) == 476
    assert len(args.moliere_tables_sha256) == 64
hydro = Path('runtime/staged_hydro/C0-5_event_00345')
assert (hydro / 'evolution_all_xyeta.dat').is_file()
out = Path(args.run_name) / 'aa/fake/task_00007'
out.mkdir(parents=True)
(out / 'HYBRID_Hadrons.out').write_text('# event 0\\nweight 1 cross 1\\nend\\n')
(out / 'summary.tsv').write_text('variant\\nno_prehydro\\n')
""".replace("PREHYDRO_ONLY", repr(prehydro_only)).replace("DO_MOLIERE", repr(do_moliere))
            )
            hydro_archive = payloads / "hydro/C0-5/event_00345.tar.gz"
            with tarfile.open(hydro_archive, "w:gz") as archive:
                base = "C0-5_event_00345"
                self.add_bytes(archive, f"{base}/evolution_all_xyeta.dat", b"hydro")
                self.add_bytes(archive, f"{base}/NcollList.dat", b"0 0\n1 1\n")
                self.add_bytes(
                    archive,
                    f"{base}/README_staged_event.txt",
                    b"hydro_slot = 12\nevent_id = 345\nncoll_positions = 2\n",
                )
            digest = hashlib.sha256(hydro_archive.read_bytes()).hexdigest()
            if not valid_checksum:
                digest = "0" * 64
            moliere_archive = moliere_payloads / "a10_tables.zip"
            with zipfile.ZipFile(moliere_archive, "w", compression=zipfile.ZIP_DEFLATED) as archive:
                for species in ("quark_tables", "gluon_tables"):
                    for index in range(476):
                        archive.writestr(
                            f"a10_tables/{species}/m{index}x_g_1_d_1_n_1.dat",
                            "1\n",
                        )
            moliere_digest = hashlib.sha256(moliere_archive.read_bytes()).hexdigest()
            manifest = runtime_root / "runtime/aa_task_manifest.tsv"
            manifest.write_text(
                "task_id\thard_seed\tmilestone_block\thydro_slot\thydro_event_id\t"
                "hydro_ncoll\thydro_dir\thydro_payload_key\thydro_payload_sha256\n"
                f"7\t900007\t1\t12\t345\t2\tC0-5_event_00345\t"
                f"hydro/C0-5/event_00345.tar.gz\t{digest}\n"
            )
            with tarfile.open(
                runtime_payloads / "mmli_runtime_alma9.tar.gz", "w:gz"
            ) as archive:
                archive.add(runtime_root / "runtime", arcname="runtime")
                archive.add(runtime_root / "bin", arcname="bin")

            job = root / "job"
            job.mkdir()
            env = os.environ.copy()
            env.update(
                {
                    "PATH": f"{fake_bin}:{env['PATH']}",
                    "MOCK_EOS": str(root / "mock_eos"),
                    "EOS_BASE": "/eos/test/campaign",
                    "PAYLOAD_EOS_BASE": "/eos/test/shared",
                    "RUNTIME_PAYLOAD_EOS_BASE": "/eos/test/runtime",
                    "KIND": "aa",
                    "SEED_OFFSET": "900000",
                    "EVENTS": "1",
                    "RUN_NAME": "runs/v2_wrapper_test",
                    "RUN_PREHYDRO_PAIR": "false",
                    "RUN_PREHYDRO_ONLY": str(prehydro_only).lower(),
                    "NO_PREHYDRO_ALPHA": "0.37",
                    "PREHYDRO_ALPHA": "0.355",
                    "BROADENING_K": "15.0",
                    "DO_MOLIERE": str(do_moliere).lower(),
                    "MOLIERE_TABLES_EOS_BASE": "/eos/test/moliere",
                    "MOLIERE_TABLES_KEY": "payloads/a10_tables.zip",
                    "MOLIERE_TABLES_SHA256": moliere_digest,
                    "AA_TASK_MANIFEST": "runtime/aa_task_manifest.tsv",
                    "TOLERATE_CHUNK_FAILURE": "false",
                }
            )
            result = subprocess.run(
                ["bash", str(RUNNER), "7"],
                cwd=job,
                env=env,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=False,
            )
            result.status_text = (job / "chunk_status.txt").read_text()  # type: ignore[attr-defined]
            return result

    def test_v2_fetches_and_records_assigned_hydro_payload(self) -> None:
        result = self.run_v2_wrapper(valid_checksum=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("hydro_slot=12", result.status_text)  # type: ignore[attr-defined]
        self.assertIn("hydro_event_id=345", result.status_text)  # type: ignore[attr-defined]
        self.assertIn("hydro_ncoll=2", result.status_text)  # type: ignore[attr-defined]

    def test_v2_rejects_hydro_payload_checksum_mismatch(self) -> None:
        result = self.run_v2_wrapper(valid_checksum=False)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("FAILED", result.stdout + result.stderr)

    def test_v2_passes_prehydro_only_mode_and_records_it(self) -> None:
        result = self.run_v2_wrapper(valid_checksum=True, prehydro_only=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("run_prehydro_pair=false", result.status_text)  # type: ignore[attr-defined]
        self.assertIn("run_prehydro_only=true", result.status_text)  # type: ignore[attr-defined]

    def test_v2_stages_and_records_moliere_tables(self) -> None:
        result = self.run_v2_wrapper(valid_checksum=True, do_moliere=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("do_moliere=true", result.status_text)  # type: ignore[attr-defined]
        self.assertIn("moliere_mode=legacy_resolved", result.status_text)  # type: ignore[attr-defined]


if __name__ == "__main__":
    unittest.main()
