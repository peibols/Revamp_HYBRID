#!/usr/bin/env python3

from __future__ import annotations

import os
from pathlib import Path
import subprocess
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


if __name__ == "__main__":
    unittest.main()
