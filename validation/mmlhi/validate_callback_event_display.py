#!/usr/bin/env python3
"""Verify that the standard Moliere Apply observer saw real hard scatterings."""

from __future__ import annotations

import argparse
from pathlib import Path

import uproot


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("event_display", type=Path)
    args = parser.parse_args()

    with uproot.open(args.event_display) as root_file:
        if "DetailedRecords" not in root_file:
            raise RuntimeError("event display has no DetailedRecords tree")
        records = root_file["DetailedRecords"]
        record_types = records["record_type"].array(library="np")
        labels = records["label"].array(library="np")

    scatterings = int((record_types == "moliere_scattering").sum())
    recoilers = int((labels == "moliere_recoiler").sum())
    holes = int((labels == "moliere_hole").sum())
    if scatterings <= 0:
        raise RuntimeError("Apply observer did not receive a hard-scattering candidate")
    if recoilers != scatterings or holes != scatterings:
        raise RuntimeError(
            "Apply observer response accounting failed: "
            f"scatterings={scatterings}, recoilers={recoilers}, holes={holes}"
        )

    print(
        "Apply observer callback contract passed: "
        f"scatterings={scatterings}, recoilers={recoilers}, holes={holes}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
