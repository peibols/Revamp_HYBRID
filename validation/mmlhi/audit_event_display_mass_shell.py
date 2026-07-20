#!/usr/bin/env python3
"""Report heavy-parton mass-shell changes in an MMLI event-display ROOT file."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import uproot


def invariant_mass2(px: float, py: float, pz: float, energy: float) -> float:
    return energy * energy - px * px - py * py - pz * pz


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root_file", type=Path)
    parser.add_argument(
        "--min-abs-delta-mass2",
        type=float,
        default=1.0e-8,
        help="minimum absolute invariant-mass-squared change to print in GeV^2",
    )
    args = parser.parse_args()

    with uproot.open(args.root_file) as root_file:
        records = root_file["DetailedRecords"].arrays(
            [
                "event_id",
                "record_id",
                "parton_index",
                "pdg_id",
                "record_type",
                "t",
                "t_end",
                "temperature",
                "px",
                "py",
                "pz",
                "e",
                "px_end",
                "py_end",
                "pz_end",
                "e_end",
            ],
            library="np",
        )

    print(
        "event_id\trecord_id\tparton_index\tpdg_id\trecord_type\t"
        "t_start\tt_end\ttemperature\tmass2_start\tmass2_end\t"
        "mass_start\tmass_end\t"
        "energy_start\tenergy_end\tdelta_energy"
    )
    heavy_rows = np.flatnonzero(np.isin(np.abs(records["pdg_id"]), (4, 5)))
    for row in heavy_rows:
        mass2_start = invariant_mass2(
            records["px"][row],
            records["py"][row],
            records["pz"][row],
            records["e"][row],
        )
        mass2_end = invariant_mass2(
            records["px_end"][row],
            records["py_end"][row],
            records["pz_end"][row],
            records["e_end"][row],
        )
        if abs(mass2_end - mass2_start) < args.min_abs_delta_mass2:
            continue
        mass_start = math.sqrt(max(0.0, mass2_start))
        mass_end = math.sqrt(max(0.0, mass2_end))
        print(
            f"{records['event_id'][row]}\t{records['record_id'][row]}\t"
            f"{records['parton_index'][row]}\t{records['pdg_id'][row]}\t"
            f"{records['record_type'][row]}\t{records['t'][row]:.12g}\t"
            f"{records['t_end'][row]:.12g}\t{records['temperature'][row]:.12g}\t"
            f"{mass2_start:.12g}\t{mass2_end:.12g}\t"
            f"{mass_start:.12g}\t{mass_end:.12g}\t{records['e'][row]:.12g}\t"
            f"{records['e_end'][row]:.12g}\t"
            f"{records['e_end'][row] - records['e'][row]:.12g}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
