#!/usr/bin/env python3
"""Produce strict three-way OO RAA, jet-spectrum, and substructure overlays."""

from __future__ import annotations

import argparse
import awkward as ak
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot

import analyze_oo_jet_raa as jet_raa
import convert_oo_paired_to_root as paired


VARIANTS = ("no_prehydro", "prehydro_alpha037", "prehydro_alpha0335")
STYLES = {
    "no_prehydro": {
        "label": "No pre-hydro",
        "color": "#0072B2",
        "marker": "o",
    },
    "prehydro_alpha037": {
        "label": r"Plan B, $\alpha=0.37$",
        "color": "#D55E00",
        "marker": "s",
    },
    "prehydro_alpha0335": {
        "label": r"Plan B, $\alpha=0.335$",
        "color": "#009E73",
        "marker": "^",
    },
}
RADII = (0.1, 0.2, 0.4, 0.8)
CHARGED_ABS = (211, 321, 2212)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        "nan"
                        if isinstance(value, float) and not math.isfinite(value)
                        else f"{value:.12e}"
                        if isinstance(value, float)
                        else value
                    )
                    for key, value in row.items()
                }
            )


def aligned_pair_metadata(
    alpha037_root: Path, alpha0335_root: Path, expected: int
) -> dict[str, np.ndarray]:
    branches = [
        "pairId",
        "chunkId",
        "hydroIndex",
        "seed",
        "eventNumber",
        "eventWeight",
        "sigmaGen",
        "hardX",
        "hardY",
    ]
    with uproot.open(alpha037_root) as root_file:
        first = root_file["Pairs"].arrays(branches, library="np")
    with uproot.open(alpha0335_root) as root_file:
        second = root_file["Pairs"].arrays(branches, library="np")
    if len(first["pairId"]) != expected:
        raise ValueError(
            f"ROOT pair count={len(first['pairId'])}, expected {expected}"
        )
    for branch in branches:
        if not np.array_equal(first[branch], second[branch]):
            raise ValueError(f"aligned ROOT files differ in Pairs/{branch}")
    if len(np.unique(first["pairId"])) != expected:
        raise ValueError("ROOT pair IDs are not unique")
    if len(np.unique(first["seed"])) != expected:
        raise ValueError("ROOT seeds are not unique")
    return first


def hadron_matrix(
    root_path: Path,
    tree_name: str,
    pair_ids: np.ndarray,
    weights: np.ndarray,
    bins: np.ndarray,
) -> np.ndarray:
    branches = [
        "pairId",
        "hadronPt",
        "hadronEta",
        "hadronStatus",
        "hadronID",
    ]
    with uproot.open(root_path) as root_file:
        arrays = root_file[tree_name].arrays(branches, library="ak")
    if not np.array_equal(ak.to_numpy(arrays.pairId), pair_ids):
        raise ValueError(f"{root_path}:{tree_name} pair ordering mismatch")
    charged = (
        (abs(arrays.hadronID) == CHARGED_ABS[0])
        | (abs(arrays.hadronID) == CHARGED_ABS[1])
        | (abs(arrays.hadronID) == CHARGED_ABS[2])
    )
    accepted = charged & (abs(arrays.hadronEta) < 1.0)
    signed = ak.where(arrays.hadronStatus < 0, -1.0, 1.0)
    matrix = np.zeros((len(pair_ids), len(bins) - 1), dtype=np.float64)
    for index, (low, high) in enumerate(zip(bins[:-1], bins[1:])):
        selected = accepted & (arrays.hadronPt >= low) & (arrays.hadronPt < high)
        counts = ak.to_numpy(ak.sum(signed[selected], axis=1)).astype(np.float64)
        matrix[:, index] = weights * counts / (high - low)
    return matrix


def load_pp_hadron(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    rows = [row for row in read_tsv(path) if row["variant"] == "no_prehydro"]
    if not rows:
        raise ValueError("pp hadron table contains no no_prehydro rows")
    bins = np.array(
        [float(rows[0]["pt_low"]), *[float(row["pt_high"]) for row in rows]],
        dtype=np.float64,
    )
    values = np.array([float(row["pp_spectrum"]) for row in rows])
    errors = np.array([float(row["pp_spectrum_stat_err"]) for row in rows])
    pp_events = {int(row["pp_events"]) for row in rows}
    if pp_events != {1_000_000}:
        raise ValueError(f"unexpected pp event counts: {sorted(pp_events)}")
    return bins, values, errors, pp_events.pop()


def make_hadron_raa(
    *,
    alpha037_root: Path,
    alpha0335_root: Path,
    pairs: dict[str, np.ndarray],
    pp_table: Path,
    out_dir: Path,
) -> tuple[list[dict[str, object]], Path]:
    bins, pp_values, pp_errors, pp_events = load_pp_hadron(pp_table)
    pair_ids = pairs["pairId"]
    weights = pairs["eventWeight"].astype(np.float64)
    sigma_gen = pairs["sigmaGen"].astype(np.float64)
    matrices = {
        "no_prehydro": hadron_matrix(
            alpha037_root, "noPrehydro/Hadrons", pair_ids, weights, bins
        ),
        "prehydro_alpha037": hadron_matrix(
            alpha037_root, "withPrehydro/Hadrons", pair_ids, weights, bins
        ),
        "prehydro_alpha0335": hadron_matrix(
            alpha0335_root, "withPrehydro/Hadrons", pair_ids, weights, bins
        ),
    }
    rows: list[dict[str, object]] = []
    for variant in VARIANTS:
        aa_values, aa_errors, weight_sum, sigma_merged = jet_raa.aggregate_spectrum(
            matrices[variant], weights, sigma_gen
        )
        raa, raa_errors = jet_raa.independent_ratio(
            aa_values, aa_errors, pp_values, pp_errors
        )
        for index, (low, high) in enumerate(zip(bins[:-1], bins[1:])):
            rows.append(
                {
                    "variant": variant,
                    "pt_low": low,
                    "pt_high": high,
                    "pt_center": math.sqrt(low * high),
                    "raa": raa[index],
                    "stat_error": raa_errors[index],
                    "aa_spectrum_mb_per_gev": aa_values[index],
                    "aa_stat_error_mb_per_gev": aa_errors[index],
                    "pp_spectrum_mb_per_gev": pp_values[index],
                    "pp_stat_error_mb_per_gev": pp_errors[index],
                    "aa_events": len(pair_ids),
                    "pp_events": pp_events,
                    "aa_weight_sum": weight_sum,
                    "aa_sigma_merged_mb": sigma_merged,
                }
            )
    table = out_dir / "oo5360_v3_hadron_raa.tsv"
    write_tsv(table, rows)

    figure, axis = plt.subplots(figsize=(7.4, 5.1))
    for variant in VARIANTS:
        selected = [row for row in rows if row["variant"] == variant]
        x = np.array([float(row["pt_center"]) for row in selected])
        low = np.array([float(row["pt_low"]) for row in selected])
        high = np.array([float(row["pt_high"]) for row in selected])
        style = STYLES[variant]
        axis.errorbar(
            x,
            [float(row["raa"]) for row in selected],
            xerr=[x - low, high - x],
            yerr=[float(row["stat_error"]) for row in selected],
            color=style["color"],
            marker=style["marker"],
            label=style["label"],
            linewidth=1.6,
            markersize=5.0,
            capsize=2.0,
        )
    axis.axhline(1.0, color="0.5", linewidth=1.0)
    axis.set_xscale("log")
    axis.set_xlim(3.8, 160)
    axis.set_ylim(0.0, 1.5)
    axis.set_xlabel(r"charged hadron $p_T$ [GeV]")
    axis.set_ylabel(r"$R_{\rm AA}$")
    axis.set_title(r"O+O 5.36 TeV, 0--5%, $|\eta|<1$, matched tasks")
    axis.grid(alpha=0.22)
    axis.legend(frameon=False, fontsize=9)
    figure.tight_layout()
    plot = out_dir / "oo5360_v3_hadron_raa.pdf"
    figure.savefig(plot)
    figure.savefig(plot.with_suffix(".png"), dpi=180)
    plt.close(figure)
    return rows, plot


def merge_jet_raa(
    alpha037_path: Path, alpha0335_path: Path, out_dir: Path
) -> tuple[list[dict[str, object]], list[Path]]:
    tables = {
        "alpha037": read_tsv(alpha037_path),
        "alpha0335": read_tsv(alpha0335_path),
    }
    key = lambda row: (float(row["radius"]), float(row["pt_low"]))
    no037 = {key(row): row for row in tables["alpha037"] if row["variant"] == "noPrehydro"}
    no335 = {key(row): row for row in tables["alpha0335"] if row["variant"] == "noPrehydro"}
    if set(no037) != set(no335):
        raise ValueError("jet RAA no-prehydro bin sets differ")
    for bin_key in no037:
        for field in (
            "raa",
            "stat_error",
            "aa_spectrum_mb_per_gev",
            "aa_stat_error_mb_per_gev",
            "aa_events",
        ):
            if not math.isclose(
                float(no037[bin_key][field]),
                float(no335[bin_key][field]),
                rel_tol=1e-12,
                abs_tol=1e-15,
            ):
                raise ValueError(f"jet RAA no-prehydro mismatch in {bin_key} {field}")

    rows: list[dict[str, object]] = []
    selections = (
        ("no_prehydro", tables["alpha037"], "noPrehydro"),
        ("prehydro_alpha037", tables["alpha037"], "withPrehydro"),
        ("prehydro_alpha0335", tables["alpha0335"], "withPrehydro"),
    )
    for output_variant, table, input_variant in selections:
        for row in table:
            if row["variant"] != input_variant:
                continue
            converted: dict[str, object] = dict(row)
            converted["variant"] = output_variant
            rows.append(converted)
    table_path = out_dir / "oo5360_v3_jet_raa.tsv"
    write_tsv(table_path, rows)

    plots: list[Path] = []
    figure, axes = plt.subplots(2, 2, figsize=(10.2, 7.4), sharex=True, sharey=True)
    for axis, radius in zip(axes.flat, RADII):
        for variant in VARIANTS:
            selected = [
                row
                for row in rows
                if row["variant"] == variant
                and math.isclose(float(row["radius"]), radius)
            ]
            selected.sort(key=lambda row: float(row["pt_low"]))
            style = STYLES[variant]
            axis.errorbar(
                [float(row["pt_center"]) for row in selected],
                [float(row["raa"]) for row in selected],
                yerr=[float(row["stat_error"]) for row in selected],
                color=style["color"],
                marker=style["marker"],
                label=style["label"],
                linewidth=1.4,
                markersize=4.2,
                capsize=1.8,
            )
        axis.axhline(1.0, color="0.55", linewidth=0.9)
        axis.set_xscale("log")
        axis.set_xlim(15.0, 230.0)
        axis.set_ylim(0.35, 1.35)
        axis.set_title(f"anti-$k_T$ R={radius:.1f}")
        axis.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False, fontsize=8)
    for axis in axes[-1, :]:
        axis.set_xlabel(r"jet $p_T$ [GeV]")
    for axis in axes[:, 0]:
        axis.set_ylabel(r"jet $R_{\rm AA}$")
    figure.suptitle(r"O+O 5.36 TeV, 0--5%, matched tasks; $|\eta_{\rm jet}|<2$")
    figure.tight_layout()
    plot = out_dir / "oo5360_v3_jet_raa.pdf"
    figure.savefig(plot)
    figure.savefig(plot.with_suffix(".png"), dpi=180)
    plt.close(figure)
    plots.append(plot)
    return rows, plots


def mapped_rows(
    alpha037_path: Path,
    alpha0335_path: Path,
) -> tuple[list[dict[str, str]], list[dict[str, str]], list[dict[str, str]]]:
    first = read_tsv(alpha037_path)
    second = read_tsv(alpha0335_path)
    output: list[dict[str, str]] = []
    for output_variant, rows, input_variant in (
        ("no_prehydro", first, "noPrehydro"),
        ("prehydro_alpha037", first, "withPrehydro"),
        ("prehydro_alpha0335", second, "withPrehydro"),
    ):
        for row in rows:
            if row.get("variant") != input_variant:
                continue
            converted = dict(row)
            converted["variant"] = output_variant
            output.append(converted)
    return output, first, second


def plot_substructure(
    *,
    hist037: Path,
    hist0335: Path,
    summary037: Path,
    summary0335: Path,
    paired_observables037: Path,
    paired_observables0335: Path,
    paired_failures037: Path,
    paired_failures0335: Path,
    out_dir: Path,
) -> tuple[list[Path], dict[str, object]]:
    hist_rows, hist_first, hist_second = mapped_rows(hist037, hist0335)
    summary_rows, _, _ = mapped_rows(summary037, summary0335)
    write_tsv(out_dir / "oo5360_v3_jet_variables_pt30_histograms.tsv", hist_rows)
    write_tsv(out_dir / "oo5360_v3_jet_variables_pt30_summary.tsv", summary_rows)

    plots: list[Path] = []
    figure, axes = plt.subplots(
        2,
        4,
        figsize=(14.4, 6.2),
        sharex="col",
        gridspec_kw={"height_ratios": [3.0, 1.0], "hspace": 0.06, "wspace": 0.28},
    )
    for column, radius in enumerate(RADII):
        upper = axes[0, column]
        lower = axes[1, column]
        selected_by_variant: dict[str, list[dict[str, str]]] = {}
        for variant in VARIANTS:
            selected = [
                row
                for row in hist_rows
                if row["variant"] == variant
                and row["variable"] == "ptd"
                and row["bin_kind"] == "physical"
                and math.isclose(float(row["radius"]), radius)
            ]
            selected.sort(key=lambda row: int(row["bin_index"]))
            selected_by_variant[variant] = selected
            widths = np.array(
                [float(row["bin_high"]) - float(row["bin_low"]) for row in selected]
            )
            values = np.array(
                [float(row["differential_cross_section_mb"]) for row in selected]
            )
            errors = np.array([float(row["stat_error_mb"]) for row in selected])
            norm = float(np.sum(values * widths))
            style = STYLES[variant]
            upper.errorbar(
                [float(row["bin_center"]) for row in selected],
                values / norm,
                yerr=errors / norm,
                color=style["color"],
                marker=style["marker"],
                label=style["label"],
                linewidth=1.3,
                markersize=3.7,
                capsize=1.5,
            )

        upper.set_title(f"R={radius:.1f}")
        upper.grid(alpha=0.2)
        upper.tick_params(labelbottom=False)
        if column == 0:
            upper.set_ylabel("normalized weighted density")
            upper.legend(frameon=False, fontsize=8)

        lower.axhline(1.0, color="0.45", linewidth=0.9)
        ratio_values: list[np.ndarray] = []
        ratio_errors: list[np.ndarray] = []
        no_rows = selected_by_variant["no_prehydro"]
        no_values = np.array(
            [float(row["differential_cross_section_mb"]) for row in no_rows]
        )
        no_errors = np.array([float(row["stat_error_mb"]) for row in no_rows])
        denominator_is_resolved = no_values > 2.0 * no_errors
        for variant in ("prehydro_alpha037", "prehydro_alpha0335"):
            selected = selected_by_variant[variant]
            ratio = np.array([float(row["pre_over_no_ratio"]) for row in selected])
            error = np.array([float(row["ratio_stat_error"]) for row in selected])
            centers = np.array([float(row["bin_center"]) for row in selected])
            finite = (
                np.isfinite(ratio)
                & np.isfinite(error)
                & denominator_is_resolved
                & (error < 0.5)
            )
            style = STYLES[variant]
            lower.errorbar(
                centers[finite],
                ratio[finite],
                yerr=error[finite],
                color=style["color"],
                marker=style["marker"],
                linestyle="none",
                markersize=3.0,
                linewidth=0.8,
                capsize=1.2,
            )
            ratio_values.append(ratio[finite])
            ratio_errors.append(error[finite])
        combined_ratio = np.concatenate(ratio_values)
        combined_error = np.concatenate(ratio_errors)
        low = min(1.0, float(np.min(combined_ratio - combined_error)))
        high = max(1.0, float(np.max(combined_ratio + combined_error)))
        padding = max(0.04, 0.12 * (high - low))
        lower.set_ylim(max(0.0, low - padding), high + padding)
        lower.set_xlabel(r"$p_{TD}$")
        lower.grid(alpha=0.2)
        if column == 0:
            lower.set_ylabel("Plan B / no")

    figure.suptitle(r"Independent jets, $p_T>30$ GeV; matched hard-event set")
    figure.text(
        0.5,
        0.012,
        "Upper: unit-normalized weighted distributions. Lower: paired differential-yield "
        "ratios to no pre-hydro with delete-one-event jackknife errors; ratio points require "
        r"no-pre-hydro bin content $>2\sigma$ and ratio error $<0.5$.",
        ha="center",
        fontsize=8.5,
    )
    figure.subplots_adjust(top=0.90, bottom=0.13, left=0.065, right=0.99)
    ptd_plot = out_dir / "oo5360_v3_ptd_distributions.pdf"
    figure.savefig(ptd_plot)
    figure.savefig(ptd_plot.with_suffix(".png"), dpi=180)
    plt.close(figure)
    plots.append(ptd_plot)

    figure, axes = plt.subplots(2, 2, figsize=(10.2, 7.4), sharex=True, sharey=True)
    for axis, radius in zip(axes.flat, RADII):
        for label, rows, color, marker in (
            (r"$\alpha=0.37$ / no pre-hydro", hist_first, "#D55E00", "s"),
            (r"$\alpha=0.335$ / no pre-hydro", hist_second, "#009E73", "^"),
        ):
            selected = [
                row
                for row in rows
                if row["variant"] == "withPrehydro"
                and row["variable"] == "pt"
                and row["bin_kind"] == "physical"
                and math.isclose(float(row["radius"]), radius)
            ]
            selected.sort(key=lambda row: int(row["bin_index"]))
            axis.errorbar(
                [float(row["bin_center"]) for row in selected],
                [float(row["pre_over_no_ratio"]) for row in selected],
                yerr=[float(row["ratio_stat_error"]) for row in selected],
                color=color,
                marker=marker,
                label=label,
                linewidth=1.3,
                markersize=3.7,
                capsize=1.5,
            )
        axis.axhline(1.0, color="0.55", linewidth=0.9)
        axis.set_xscale("log")
        axis.set_xlim(28.0, 300.0)
        axis.set_ylim(0.78, 1.12)
        axis.set_title(f"R={radius:.1f}")
        axis.grid(alpha=0.2)
    axes[0, 0].legend(frameon=False, fontsize=8)
    for axis in axes[-1, :]:
        axis.set_xlabel(r"jet $p_T$ [GeV]")
    for axis in axes[:, 0]:
        axis.set_ylabel("spectrum ratio to no pre-hydro")
    figure.suptitle(r"Same-seed paired-weight jet spectra, $30<p_T<300$ GeV")
    figure.tight_layout()
    spectrum_plot = out_dir / "oo5360_v3_jet_spectrum_ratios.pdf"
    figure.savefig(spectrum_plot)
    figure.savefig(spectrum_plot.with_suffix(".png"), dpi=180)
    plt.close(figure)
    plots.append(spectrum_plot)

    def selected_observable(path: Path, observable: str) -> dict[float, dict[str, str]]:
        rows = read_tsv(path)
        return {
            float(row["radius"]): row
            for row in rows
            if row["flavor"] == "all"
            and row["category"] == "all"
            and row["observable"] == observable
        }

    def selected_failure(path: Path) -> dict[float, dict[str, str]]:
        rows = read_tsv(path)
        return {
            float(row["radius"]): row
            for row in rows
            if row["flavor"] == "all"
            and row["observable"] == "delta_fail_fraction"
        }

    shifts = {
        "prehydro_alpha037": {
            "ptd": selected_observable(
                paired_observables037, "relative_mean_ptd_shift"
            ),
            "fail": selected_failure(paired_failures037),
        },
        "prehydro_alpha0335": {
            "ptd": selected_observable(
                paired_observables0335, "relative_mean_ptd_shift"
            ),
            "fail": selected_failure(paired_failures0335),
        },
    }
    figure, axes = plt.subplots(1, 2, figsize=(10.0, 4.3))
    for variant in ("prehydro_alpha037", "prehydro_alpha0335"):
        style = STYLES[variant]
        axes[0].errorbar(
            RADII,
            [100.0 * float(shifts[variant]["ptd"][radius]["value"]) for radius in RADII],
            yerr=[
                100.0
                * float(shifts[variant]["ptd"][radius]["eventJackknifeError"])
                for radius in RADII
            ],
            color=style["color"],
            marker=style["marker"],
            label=style["label"],
            linewidth=1.5,
            capsize=2.2,
        )
        axes[1].errorbar(
            RADII,
            [100.0 * float(shifts[variant]["fail"][radius]["value"]) for radius in RADII],
            yerr=[
                100.0
                * float(shifts[variant]["fail"][radius]["eventJackknifeError"])
                for radius in RADII
            ],
            color=style["color"],
            marker=style["marker"],
            label=style["label"],
            linewidth=1.5,
            capsize=2.2,
        )
    for axis in axes:
        axis.axhline(0.0, color="0.55", linewidth=0.9)
        axis.set_xticks(RADII)
        axis.set_xlabel("jet radius R")
        axis.grid(alpha=0.2)
    axes[0].set_ylabel(r"matched mean-$p_{TD}$ shift [\%]")
    axes[1].set_ylabel("matched SD-fail change [percentage points]")
    axes[0].legend(frameon=False, fontsize=8)
    figure.suptitle(r"No-prehydro-selected matched jets, $p_T>30$ GeV, $|\eta|<2$")
    figure.tight_layout()
    shift_plot = out_dir / "oo5360_v3_matched_substructure_shifts.pdf"
    figure.savefig(shift_plot)
    figure.savefig(shift_plot.with_suffix(".png"), dpi=180)
    plt.close(figure)
    plots.append(shift_plot)

    summary = {
        variant: {
            f"R{radius:.1f}": {
                "relativeMeanPtdShift": float(shifts[variant]["ptd"][radius]["value"]),
                "relativeMeanPtdShiftError": float(
                    shifts[variant]["ptd"][radius]["eventJackknifeError"]
                ),
                "deltaSoftDropFailFraction": float(
                    shifts[variant]["fail"][radius]["value"]
                ),
                "deltaSoftDropFailFractionError": float(
                    shifts[variant]["fail"][radius]["eventJackknifeError"]
                ),
            }
            for radius in RADII
        }
        for variant in ("prehydro_alpha037", "prehydro_alpha0335")
    }
    return plots, summary


def row_at(
    rows: list[dict[str, object]], variant: str, low: float, radius: float | None = None
) -> dict[str, object]:
    selected = [
        row
        for row in rows
        if row["variant"] == variant
        and math.isclose(float(row["pt_low"]), low)
        and (radius is None or math.isclose(float(row["radius"]), radius))
    ]
    if len(selected) != 1:
        raise ValueError(f"expected one row for {variant}, low={low}, radius={radius}")
    return selected[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root-alpha037", type=Path, required=True)
    parser.add_argument("--root-alpha0335", type=Path, required=True)
    parser.add_argument("--pp-hadron-table", type=Path, required=True)
    parser.add_argument("--jet-raa-alpha037", type=Path, required=True)
    parser.add_argument("--jet-raa-alpha0335", type=Path, required=True)
    parser.add_argument("--variables-hist-alpha037", type=Path, required=True)
    parser.add_argument("--variables-hist-alpha0335", type=Path, required=True)
    parser.add_argument("--variables-summary-alpha037", type=Path, required=True)
    parser.add_argument("--variables-summary-alpha0335", type=Path, required=True)
    parser.add_argument("--paired-observables-alpha037", type=Path, required=True)
    parser.add_argument("--paired-observables-alpha0335", type=Path, required=True)
    parser.add_argument("--paired-failures-alpha037", type=Path, required=True)
    parser.add_argument("--paired-failures-alpha0335", type=Path, required=True)
    parser.add_argument("--expected-events", type=int, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    pairs = aligned_pair_metadata(
        args.root_alpha037, args.root_alpha0335, args.expected_events
    )
    hadron_rows, hadron_plot = make_hadron_raa(
        alpha037_root=args.root_alpha037,
        alpha0335_root=args.root_alpha0335,
        pairs=pairs,
        pp_table=args.pp_hadron_table,
        out_dir=args.out_dir,
    )
    jet_rows, jet_plots = merge_jet_raa(
        args.jet_raa_alpha037, args.jet_raa_alpha0335, args.out_dir
    )
    substructure_plots, substructure = plot_substructure(
        hist037=args.variables_hist_alpha037,
        hist0335=args.variables_hist_alpha0335,
        summary037=args.variables_summary_alpha037,
        summary0335=args.variables_summary_alpha0335,
        paired_observables037=args.paired_observables_alpha037,
        paired_observables0335=args.paired_observables_alpha0335,
        paired_failures037=args.paired_failures_alpha037,
        paired_failures0335=args.paired_failures_alpha0335,
        out_dir=args.out_dir,
    )

    key_results = {
        "status": "PASS",
        "publicationStatus": "PROVISIONAL_COMPLETION_ORDER_MATCHED_DIAGNOSTIC",
        "matchedTriplets": args.expected_events,
        "normalization": "PythiaParallel sigmaGen/sum(weight), applied once; common 1M pp denominator",
        "uncertainty": "delete-one-AA-event jackknife; pp uncertainty combined independently for RAA",
        "ptdRatioDisplay": {
            "denominator": "no_prehydro",
            "quantity": "paired differential-yield ratio",
            "minimumNoPrehydroSignificance": 2.0,
            "maximumDisplayedRatioStatError": 0.5,
            "note": "display cuts affect ratio markers only; tables and upper distributions are unfiltered",
        },
        "hadronRaa10to14": {
            variant: {
                "value": float(row_at(hadron_rows, variant, 10.0)["raa"]),
                "error": float(row_at(hadron_rows, variant, 10.0)["stat_error"]),
            }
            for variant in VARIANTS
        },
        "jetRaaR01_30to40": {
            variant: {
                "value": float(row_at(jet_rows, variant, 30.0, 0.1)["raa"]),
                "error": float(row_at(jet_rows, variant, 30.0, 0.1)["stat_error"]),
            }
            for variant in VARIANTS
        },
        "matchedSubstructure": substructure,
        "inputs": {
            "rootAlpha037": str(args.root_alpha037.resolve()),
            "rootAlpha037Sha256": paired.sha256(args.root_alpha037),
            "rootAlpha0335": str(args.root_alpha0335.resolve()),
            "rootAlpha0335Sha256": paired.sha256(args.root_alpha0335),
        },
        "plots": [str(hadron_plot), *map(str, jet_plots), *map(str, substructure_plots)],
    }
    metadata = args.out_dir / "oo5360_v3_comparison_metadata.json"
    metadata.write_text(json.dumps(key_results, indent=2, sort_keys=True) + "\n")
    print(json.dumps(key_results, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
