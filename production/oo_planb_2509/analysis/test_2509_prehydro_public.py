#!/usr/bin/env python3
"""Regression tests for the public arXiv:2509.19430 prehydro reproduction."""

from __future__ import annotations

import csv
import hashlib
from importlib.util import module_from_spec, spec_from_file_location
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


WORK = Path(__file__).resolve().parents[1]
RUNNER_PATH = WORK / "cern_support/run_oo_validation_chunk.py"
ATTRACTOR_PATH = (
    WORK / "reference_data/qcd_kinetic_attractor_lambda10_Cinf0p87.tsv"
)
HYDRO_PATH = WORK / "staged_hydro/C0-5_idx0/evolution_all_xyeta.dat"

SPEC = spec_from_file_location("oo_prehydro_runner", RUNNER_PATH)
assert SPEC is not None and SPEC.loader is not None
RUNNER = module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


class PublicAttractorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.label, attractor = RUNNER.load_attractor_table(ATTRACTOR_PATH)
        cls.attractor = staticmethod(attractor)

    def test_published_table_checksum_and_values(self) -> None:
        digest = hashlib.sha256(ATTRACTOR_PATH.read_bytes()).hexdigest()
        self.assertEqual(digest, RUNNER.PUBLISHED_QCD_ATTRACTOR_SHA256)
        self.assertIn(digest, self.label)
        expected = {
            0.20: 0.534722025,
            0.30: 0.611917251,
            0.50: 0.707616327,
        }
        for omega, value in expected.items():
            self.assertAlmostEqual(self.attractor(omega), value, places=8)

    def test_table_rejects_out_of_range_and_modified_data(self) -> None:
        with self.assertRaises(ValueError):
            self.attractor(0.01)
        with self.assertRaises(ValueError):
            self.attractor(10.0)
        with tempfile.TemporaryDirectory() as tmp:
            modified = Path(tmp) / "modified.tsv"
            modified.write_bytes(ATTRACTOR_PATH.read_bytes() + b"\n")
            with self.assertRaisesRegex(ValueError, "SHA256"):
                RUNNER.load_attractor_table(modified)

    def test_temperature_satisfies_published_implicit_equation(self) -> None:
        tau = 0.24
        tau_hyd = 0.4
        temperature_hyd = 0.25
        eta_over_s = 0.12
        temperature, omega, e_value = RUNNER.attractor_temperature(
            tau=tau,
            tau_hyd=tau_hyd,
            temperature_hyd=temperature_hyd,
            eta_over_s=eta_over_s,
            attractor=self.attractor,
            viscous_anchor=True,
        )
        tau_natural = tau / RUNNER.HBARC_GEV_FM
        tau_hyd_natural = tau_hyd / RUNNER.HBARC_GEV_FM
        lhs = (tau_natural ** (1.0 / 3.0) * temperature) ** 4
        anchor_temperature = temperature_hyd + 2.0 * eta_over_s / (
            3.0 * tau_hyd_natural
        )
        rhs = (
            tau_hyd_natural ** (1.0 / 3.0) * anchor_temperature
        ) ** 4 * e_value
        self.assertAlmostEqual(lhs / rhs, 1.0, places=9)
        self.assertAlmostEqual(
            omega,
            tau_natural * temperature / (4.0 * math.pi * eta_over_s),
            places=12,
        )

    def test_dense_grid_linear_interpolation_error_is_below_one_per_ten_thousand(
        self,
    ) -> None:
        worst = 0.0
        for temperature_hyd in (0.145, 0.2, 0.3, 0.5, 0.8):
            for index in range(24, 39):
                tau_low = index / 100.0
                tau_high = (index + 1) / 100.0
                tau_mid = 0.5 * (tau_low + tau_high)
                low = RUNNER.attractor_temperature(
                    tau=tau_low,
                    tau_hyd=0.4,
                    temperature_hyd=temperature_hyd,
                    eta_over_s=0.12,
                    attractor=self.attractor,
                    viscous_anchor=True,
                )[0]
                high = RUNNER.attractor_temperature(
                    tau=tau_high,
                    tau_hyd=0.4,
                    temperature_hyd=temperature_hyd,
                    eta_over_s=0.12,
                    attractor=self.attractor,
                    viscous_anchor=True,
                )[0]
                direct = RUNNER.attractor_temperature(
                    tau=tau_mid,
                    tau_hyd=0.4,
                    temperature_hyd=temperature_hyd,
                    eta_over_s=0.12,
                    attractor=self.attractor,
                    viscous_anchor=True,
                )[0]
                worst = max(worst, abs((low + high) / (2.0 * direct) - 1.0))
        self.assertLess(worst, 1e-4)

    def test_generated_table_applies_onset_flow_and_conformal_eos(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "prehydro.tsv"
            RUNNER.generate_planb_prehydro(
                reference_hydro=HYDRO_PATH,
                output_path=output,
                event_id=9525,
                tau_grid=[0.20, 0.24, 0.25, 0.399, 0.40],
                tau_min=0.24,
                eta_over_s=0.12,
                eos_factor=RUNNER.QCD_CONFORMAL_EOS_FACTOR,
                attractor_label=self.label,
                attractor=self.attractor,
                viscous_anchor=True,
            )
            text = output.read_text()
            self.assertIn(RUNNER.PUBLISHED_QCD_ATTRACTOR_SHA256, text)
            self.assertIn("arXiv:2509.19430v2 Eqs. (2)-(4)", text)
            self.assertIn(RUNNER.PUBLISHED_QCD_ATTRACTOR_DOI, text)
            self.assertIn("# eos_effective_degrees = 47.5", text)
            rows = list(
                csv.reader(
                    (line for line in text.splitlines() if not line.startswith("#")),
                    delimiter="\t",
                )
            )
            self.assertEqual({float(row[1]) for row in rows}, {0.24, 0.25, 0.399})
            self.assertTrue(all(int(row[10]) == 1 for row in rows))

            row = rows[0]
            tau = float(row[1])
            temperature = float(row[6])
            energy_density = float(row[7])
            vx = float(row[8])
            vy = float(row[9])
            vx_hyd = float(row[12])
            vy_hyd = float(row[13])
            tau_hyd = RUNNER.parse_hydro_header(HYDRO_PATH)["tau0_fm"]
            self.assertAlmostEqual(
                energy_density,
                RUNNER.QCD_CONFORMAL_EOS_FACTOR * temperature**4,
                places=8,
            )
            self.assertAlmostEqual(vx, vx_hyd * tau / tau_hyd, places=7)
            self.assertAlmostEqual(vy, vy_hyd * tau / tau_hyd, places=7)

    def test_paired_cli_requires_the_published_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result = subprocess.run(
                [
                    sys.executable,
                    str(RUNNER_PATH),
                    "--kind", "aa",
                    "--task-id", "0",
                    "--seed-offset", "1",
                    "--events", "1",
                    "--run-name", "test",
                    "--run-prehydro-pair",
                ],
                cwd=tmp,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("--prehydro-attractor-table is required", result.stdout)

    def test_paired_cli_pins_published_eta_over_s_and_viscous_anchor(self) -> None:
        common = [
            sys.executable,
            str(RUNNER_PATH),
            "--kind", "aa",
            "--task-id", "0",
            "--seed-offset", "1",
            "--events", "1",
            "--run-name", "test",
            "--run-prehydro-pair",
            "--prehydro-attractor-table", str(ATTRACTOR_PATH),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            wrong_eta = subprocess.run(
                [*common, "--prehydro-eta-over-s", "0.16"],
                cwd=tmp,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(wrong_eta.returncode, 2)
            self.assertIn("must be 0.12", wrong_eta.stdout)

            no_anchor = subprocess.run(
                [*common, "--no-prehydro-viscous-anchor"],
                cwd=tmp,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                check=False,
            )
            self.assertEqual(no_anchor.returncode, 2)
            self.assertIn("incompatible", no_anchor.stdout)


if __name__ == "__main__":
    unittest.main()
