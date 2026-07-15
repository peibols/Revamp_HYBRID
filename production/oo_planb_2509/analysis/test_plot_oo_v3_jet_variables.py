#!/usr/bin/env python3
"""Regression tests for three-way V3 jet-variable plotting."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


MODULE_PATH = Path(__file__).with_name("plot_oo_v3_jet_variables.py")
sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("plot_oo_v3_jet_variables", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
v3 = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = v3
SPEC.loader.exec_module(v3)
base = v3.base


class V3JetVariableRebinTest(unittest.TestCase):
    @staticmethod
    def rows(spec: base.VariableSpec) -> list[dict[str, str]]:
        return [
            {
                "variant": "no_prehydro",
                "variable": spec.key,
                "radius": "0.4",
                "bin_index": str(index),
                "bin_low": f"{low:.12e}",
                "bin_high": f"{high:.12e}",
            }
            for index, (low, high) in enumerate(zip(spec.edges[:-1], spec.edges[1:]))
        ]

    def test_factor_two_rows_match_factor_two_schema(self) -> None:
        specs = {
            spec.key: spec
            for spec in base.variable_specs(
                4, 30.0, substructure_rebin_factor=2
            )
        }
        for key in ("zg", "rg", "mult", "totalmult", "ptd", "girth", "maxkt"):
            selected = v3.rows_for_spec(
                self.rows(specs[key]), 0.4, specs[key], "no_prehydro"
            )
            self.assertEqual(len(selected), len(specs[key].edges) - 1)

    def test_factor_two_rows_fail_unrebinned_schema(self) -> None:
        original = {spec.key: spec for spec in base.variable_specs(4, 30.0)}["ptd"]
        rebinned = {
            spec.key: spec
            for spec in base.variable_specs(
                4, 30.0, substructure_rebin_factor=2
            )
        }["ptd"]
        with self.assertRaisesRegex(ValueError, "has 10 bins, expected 20"):
            v3.rows_for_spec(
                self.rows(rebinned), 0.4, original, "no_prehydro"
            )


if __name__ == "__main__":
    unittest.main()
