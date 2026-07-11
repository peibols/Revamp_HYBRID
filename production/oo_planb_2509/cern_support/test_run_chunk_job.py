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

    def run_v2_wrapper(self, valid_checksum: bool) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            eos = root / "mock_eos/eos/test/campaign"
            payloads = eos / "payloads"
            (payloads / "hydro/C0-5").mkdir(parents=True)
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
args, _ = p.parse_known_args()
assert Path(args.aa_task_manifest).is_file()
hydro = Path('runtime/staged_hydro/C0-5_event_00345')
assert (hydro / 'evolution_all_xyeta.dat').is_file()
out = Path(args.run_name) / 'aa/fake/task_00007'
out.mkdir(parents=True)
(out / 'HYBRID_Hadrons.out').write_text('# event 0\\nweight 1 cross 1\\nend\\n')
(out / 'summary.tsv').write_text('variant\\nno_prehydro\\n')
"""
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
            manifest = runtime_root / "runtime/aa_task_manifest.tsv"
            manifest.write_text(
                "task_id\thard_seed\tmilestone_block\thydro_slot\thydro_event_id\t"
                "hydro_ncoll\thydro_dir\thydro_payload_key\thydro_payload_sha256\n"
                f"7\t900007\t1\t12\t345\t2\tC0-5_event_00345\t"
                f"hydro/C0-5/event_00345.tar.gz\t{digest}\n"
            )
            with tarfile.open(payloads / "mmli_runtime_alma9.tar.gz", "w:gz") as archive:
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
                    "KIND": "aa",
                    "SEED_OFFSET": "900000",
                    "EVENTS": "1",
                    "RUN_NAME": "runs/v2_wrapper_test",
                    "RUN_PREHYDRO_PAIR": "false",
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


if __name__ == "__main__":
    unittest.main()
