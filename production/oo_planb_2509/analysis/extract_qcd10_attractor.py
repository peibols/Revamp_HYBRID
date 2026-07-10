#!/usr/bin/env python3
"""Extract the published QCD lambda=10 energy-attractor curve."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit


N_GLUON = 16.0
N_FERMION = 36.0
ETA_OVER_S = 1.0


def fit_function(omega: np.ndarray, intercept: float, slope: float) -> np.ndarray:
    return intercept + slope / omega


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    gluon = np.loadtxt(
        args.data_root / "prl_gf_l10_xi10_exp_gluon_Tmunu_vs_time.out"
    )
    fermion = np.loadtxt(
        args.data_root / "prl_gf_l10_xi10_exp_fermion_Tmunu_vs_time.out"
    )
    if gluon.shape != fermion.shape or not np.allclose(
        gluon[:, 0], fermion[:, 0], rtol=0.0, atol=1e-12
    ):
        raise RuntimeError("gluon and fermion evolution grids do not match")

    time = gluon[:, 0]
    energy = N_GLUON * gluon[:, 8] + N_FERMION * fermion[:, 8]
    effective_degrees = N_GLUON + 7.0 * N_FERMION / 8.0
    temperature = (
        energy / np.pi**2 * 30.0 / effective_degrees
    ) ** 0.25
    omega = time * temperature / (4.0 * np.pi * ETA_OVER_S)

    fit_values = (
        effective_degrees
        * np.pi**2
        / 30.0
        * (
            time ** (1.0 / 3.0)
            * (temperature + 2.0 * ETA_OVER_S / (3.0 * time))
        )
        ** 4
    )
    fit_mask = np.isfinite(omega) & np.isfinite(fit_values) & (omega > 3.0)
    parameters, _ = curve_fit(
        fit_function, omega[fit_mask], fit_values[fit_mask]
    )
    normalization = float(parameters[0])
    attractor = energy * time ** (4.0 / 3.0) / normalization

    valid = (
        np.isfinite(omega)
        & np.isfinite(attractor)
        & (omega > 0.0)
        & (attractor > 0.0)
    )
    omega = omega[valid]
    attractor = attractor[valid]
    order = np.argsort(omega)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="ascii") as stream:
        stream.write(
            "# QCD kinetic-theory energy attractor from "
            "DOI:10.4119/unibi/2939684\n"
        )
        stream.write(
            "# Reproduces the QCD lambda=10 curve in data/figure1.py "
            "(C_infinity=0.87 label)\n"
        )
        stream.write(
            "# Source paper: arXiv:1908.02866; cited by arXiv:2509.19430\n"
        )
        stream.write(
            f"# late_time_fit_normalization = {normalization:.16g}\n"
        )
        stream.write("# omega\tE_omega\n")
        for x_value, e_value in zip(omega[order], attractor[order]):
            stream.write(f"{x_value:.12g}\t{e_value:.12g}\n")

    print(f"wrote {args.output}")
    print(
        f"points={len(omega)} omega_min={omega.min():.8g} "
        f"omega_max={omega.max():.8g}"
    )
    print(f"late_time_fit_normalization={normalization:.12g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
