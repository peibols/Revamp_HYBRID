#!/usr/bin/env python3
"""Merge two aligned pair analyses into the OO V3 formation-time comparison."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
import numpy as np

import plot_oo_jet_formation_time as pair
import plot_oo_jet_variables as base


VARIANTS = ("no_prehydro", "prehydro_alpha037", "prehydro_alpha0335")
STYLES = {
    "no_prehydro": {"label": "No pre-hydro", "color": "#111111", "marker": None},
    "prehydro_alpha037": {
        "label": r"Pre-hydro, $\alpha=0.37$",
        "color": "#0072B2",
        "marker": "o",
    },
    "prehydro_alpha0335": {
        "label": r"Pre-hydro, $\alpha=0.335$",
        "color": "#D62728",
        "marker": "o",
    },
}
SPECTRUM_STRUCTURE = (
    "radius",
    "pt_low_exclusive_GeV",
    "pt_high_inclusive_GeV",
    "splitting_selection",
    "bin_index",
    "log10_tau_low",
    "log10_tau_high",
    "log10_tau_center",
)
CORRELATION_STRUCTURE = (
    "radius",
    "pt_low_exclusive_GeV",
    "pt_high_inclusive_GeV",
    "correlation",
    "x_bin_index",
    "x_low",
    "x_high",
    "log10_tau_bin_index",
    "log10_tau_low",
    "log10_tau_high",
)
SUMMARY_STRUCTURE = (
    "radius",
    "pt_low_exclusive_GeV",
    "pt_high_inclusive_GeV",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty TSV: {path}")
    return rows


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty TSV: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def row_key(row: dict[str, str], fields: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(row[field] for field in fields)


def variant_map(
    rows: list[dict[str, str]], variant: str, fields: tuple[str, ...]
) -> dict[tuple[str, ...], dict[str, str]]:
    selected = {row_key(row, fields): row for row in rows if row["variant"] == variant}
    if not selected:
        raise ValueError(f"no rows for {variant}")
    return selected


def numeric_match(first: str, second: str) -> bool:
    try:
        first_value = float(first)
        second_value = float(second)
    except ValueError:
        return first == second
    if math.isnan(first_value) and math.isnan(second_value):
        return True
    return math.isclose(first_value, second_value, rel_tol=1.0e-12, abs_tol=1.0e-15)


def assert_rows_match(
    first: dict[str, str], second: dict[str, str], fields: tuple[str, ...], label: str
) -> None:
    for field in fields:
        if not numeric_match(first[field], second[field]):
            raise ValueError(
                f"{label} differs in {field}: {first[field]} != {second[field]}"
            )


def merge_table(
    first_rows: list[dict[str, str]],
    second_rows: list[dict[str, str]],
    structure: tuple[str, ...],
    baseline_fields: tuple[str, ...],
) -> list[dict[str, str]]:
    maps = {
        "no037": variant_map(first_rows, "noPrehydro", structure),
        "pre037": variant_map(first_rows, "withPrehydro", structure),
        "no335": variant_map(second_rows, "noPrehydro", structure),
        "pre335": variant_map(second_rows, "withPrehydro", structure),
    }
    keys = set(maps["no037"])
    if any(set(item) != keys for item in maps.values()):
        raise ValueError("aligned pair tables have different row keys")
    output: list[dict[str, str]] = []
    for key in sorted(keys):
        assert_rows_match(
            maps["no037"][key],
            maps["no335"][key],
            (*structure, *baseline_fields),
            f"no-prehydro {key}",
        )
        for variant, source in (
            ("no_prehydro", maps["no037"][key]),
            ("prehydro_alpha037", maps["pre037"][key]),
            ("prehydro_alpha0335", maps["pre335"][key]),
        ):
            converted = dict(source)
            converted["variant"] = variant
            output.append(converted)
    return output


def selection_rows(
    rows: list[dict[str, str]], radius: float, pt_low: float, pt_high: float | None
) -> list[dict[str, str]]:
    high = "inf" if pt_high is None else f"{pt_high:g}"
    return [
        row
        for row in rows
        if math.isclose(float(row["radius"]), radius)
        and math.isclose(float(row["pt_low_exclusive_GeV"]), pt_low)
        and row["pt_high_inclusive_GeV"] == high
    ]


def finite_errorbar(axis, x, y, yerr, **kwargs) -> None:
    x = np.asarray(x)
    y = np.asarray(y)
    yerr = np.asarray(yerr)
    finite = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr)
    if np.any(finite):
        axis.errorbar(x[finite], y[finite], yerr=yerr[finite], **kwargs)


def plot_spectra(
    rows: list[dict[str, str]],
    radius: float,
    pt_low: float,
    pt_high: float | None,
    output_base: Path,
    sample_label: str,
) -> None:
    selected = selection_rows(rows, radius, pt_low, pt_high)
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(11.8, 7.6),
        sharex="col",
        gridspec_kw={"height_ratios": [3.0, 1.15]},
    )
    for column, split_kind in enumerate(("all", "hardest_kt")):
        split_rows = [row for row in selected if row["splitting_selection"] == split_kind]
        no_rows = sorted(
            (row for row in split_rows if row["variant"] == "no_prehydro"),
            key=lambda row: int(row["bin_index"]),
        )
        no_raw_by_bin = np.asarray(
            [int(row["raw_entries_in_bin"]) for row in no_rows], dtype=np.int64
        )
        ratios_for_limits = []
        errors_for_limits = []
        for variant in VARIANTS:
            variant_rows = sorted(
                (row for row in split_rows if row["variant"] == variant),
                key=lambda row: int(row["bin_index"]),
            )
            x = np.asarray([float(row["log10_tau_center"]) for row in variant_rows])
            values = np.asarray(
                [float(row["selected_jet_normalized_density"]) for row in variant_rows]
            )
            errors = np.asarray(
                [
                    float(row["selected_jet_normalized_density_stat_error"])
                    for row in variant_rows
                ]
            )
            if variant == "no_prehydro":
                edges = [float(variant_rows[0]["log10_tau_low"])] + [
                    float(row["log10_tau_high"]) for row in variant_rows
                ]
                axes[0, column].stairs(
                    values,
                    edges,
                    color=STYLES[variant]["color"],
                    linewidth=1.7,
                    label=STYLES[variant]["label"],
                )
            else:
                finite_errorbar(
                    axes[0, column],
                    x,
                    values,
                    errors,
                    color=STYLES[variant]["color"],
                    marker=STYLES[variant]["marker"],
                    markersize=3.2,
                    linestyle="none",
                    linewidth=1.0,
                    capsize=1.8,
                    label=STYLES[variant]["label"],
                )
            if variant != "no_prehydro":
                ratio = np.asarray([float(row["pre_over_no_ratio"]) for row in variant_rows])
                ratio_error = np.asarray(
                    [float(row["ratio_stat_error"]) for row in variant_rows]
                )
                ratio_for_plot = np.where(no_raw_by_bin >= 20, ratio, np.nan)
                error_for_plot = np.where(no_raw_by_bin >= 20, ratio_error, np.nan)
                ratios_for_limits.append(ratio_for_plot)
                errors_for_limits.append(error_for_plot)
                finite_errorbar(
                    axes[1, column],
                    x,
                    ratio_for_plot,
                    error_for_plot,
                    color=STYLES[variant]["color"],
                    marker=STYLES[variant]["marker"],
                    markersize=3.2,
                    linewidth=1.0,
                    capsize=1.8,
                    label=STYLES[variant]["label"] + " / no",
                )
        axes[0, column].set_yscale("log")
        axes[0, column].set_ylabel(
            r"$(1/\sigma_{\rm jet})\,d\sigma_{\rm split}/d\log_{10}\tau_{\rm f}$"
        )
        axes[1, column].set_ylabel("pre-hydro / no pre-hydro")
        axes[1, column].set_xlabel(r"$\log_{10}(\tau_{\rm f}/[\mathrm{fm}/c])$")
        axes[1, column].axhline(1.0, color="0.5", linewidth=0.9)
        if ratios_for_limits:
            axes[1, column].set_ylim(
                base.ratio_limits(
                    np.concatenate(ratios_for_limits), np.concatenate(errors_for_limits)
                )
            )
        axes[0, column].set_title(
            "All C/A declusterings"
            if split_kind == "all"
            else r"global hardest-$k_T$ split per jet"
        )
        for row in range(2):
            axes[row, column].grid(alpha=0.2)
            axes[row, column].axvline(math.log10(0.1), color="#009E73", linestyle="--", linewidth=0.9)
            axes[row, column].axvline(math.log10(0.24), color="#CC79A7", linestyle=":", linewidth=1.0)
    axes[0, 0].legend(frameon=False, fontsize=9)
    axes[1, 0].legend(frameon=False, fontsize=8)
    figure.suptitle(
        rf"O+O 5.36 TeV, 0--5%, anti-$k_T$ R={radius:g}, {pair.pt_label(pt_low, pt_high)}"
        + f"\n{sample_label}; exact 3D opening angle; "
        + r"$\tau_{\rm f}=\hbar cE/Q^2$",
        fontsize=11,
    )
    figure.text(
        0.5,
        0.008,
        r"Vertical lines: $\tau=0.1$ and $0.24$ fm/$c$. Jackknife ratios require $N_{\rm no}\geq20$ per bin.",
        ha="center",
        fontsize=8,
    )
    figure.tight_layout(rect=(0.0, 0.035, 1.0, 0.93))
    for suffix in ("pdf", "png"):
        figure.savefig(output_base.with_suffix(f".{suffix}"), dpi=180)
    plt.close(figure)


def grid_from_rows(
    rows: list[dict[str, str]], variant: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    selected = [row for row in rows if row["variant"] == variant]
    x_bins = max(int(row["x_bin_index"]) for row in selected) + 1
    y_bins = max(int(row["log10_tau_bin_index"]) for row in selected) + 1
    density = np.zeros((x_bins, y_bins), dtype=float)
    raw = np.zeros((x_bins, y_bins), dtype=float)
    x_edges = np.empty(x_bins + 1, dtype=float)
    y_edges = np.empty(y_bins + 1, dtype=float)
    for row in selected:
        x_index = int(row["x_bin_index"])
        y_index = int(row["log10_tau_bin_index"])
        density[x_index, y_index] = float(row["selected_jet_normalized_density"])
        raw[x_index, y_index] = float(row["raw_entries"])
        x_edges[x_index] = float(row["x_low"])
        x_edges[x_index + 1] = float(row["x_high"])
        y_edges[y_index] = float(row["log10_tau_low"])
        y_edges[y_index + 1] = float(row["log10_tau_high"])
    return density, raw, x_edges, y_edges


def plot_correlations(
    rows: list[dict[str, str]],
    radius: float,
    pt_low: float,
    pt_high: float | None,
    output_base: Path,
    sample_label: str,
) -> None:
    selected = selection_rows(rows, radius, pt_low, pt_high)
    correlation_labels = {
        "z": r"$z=\min(E_1,E_2)/E_{\rm parent}$",
        "deltaR": r"$\log_{10}\Delta R_{12}$",
        "kt": r"$\log_{10}(k_T/\mathrm{GeV})$",
    }
    figure, axes = plt.subplots(3, 3, figsize=(13.6, 10.0), sharey=True)
    axes[0, 0].set_title(r"No pre-hydro $(1/\sigma_{\rm jet})$ density")
    axes[0, 1].set_title(r"Pre-hydro $\alpha=0.37$ / no")
    axes[0, 2].set_title(r"Pre-hydro $\alpha=0.335$ / no")
    for row_index, name in enumerate(("z", "deltaR", "kt")):
        correlation_rows = [row for row in selected if row["correlation"] == name]
        no_density, no_raw, x_edges, tau_edges = grid_from_rows(
            correlation_rows, "no_prehydro"
        )
        positive = no_density[no_density > 0.0]
        vmin = max(float(np.percentile(positive, 3.0)), float(np.min(positive)))
        vmax = float(np.max(positive))
        if vmin >= vmax:
            vmin = vmax / 10.0
        no_mesh = axes[row_index, 0].pcolormesh(
            x_edges,
            tau_edges,
            np.ma.masked_less_equal(no_density.T, 0.0),
            shading="auto",
            cmap="viridis",
            norm=LogNorm(vmin=vmin, vmax=vmax),
        )
        figure.colorbar(
            no_mesh,
            ax=axes[row_index, 0],
            pad=0.01,
            label=r"$(1/\sigma_{\rm jet})d^2\sigma_{\rm split}/dx\,d\log\tau_f$",
        )
        for column, variant in enumerate(
            ("prehydro_alpha037", "prehydro_alpha0335"), start=1
        ):
            pre_density, _, pre_x_edges, pre_tau_edges = grid_from_rows(
                correlation_rows, variant
            )
            if not np.array_equal(x_edges, pre_x_edges) or not np.array_equal(
                tau_edges, pre_tau_edges
            ):
                raise ValueError("V3 correlation binning differs")
            ratio = np.divide(
                pre_density,
                no_density,
                out=np.full_like(pre_density, np.nan),
                where=no_density > 0.0,
            )
            ratio = np.ma.masked_where(no_raw < 10, ratio)
            finite = ratio.compressed()
            distance = (
                max(0.2, min(1.0, float(np.percentile(np.abs(finite - 1.0), 95.0))))
                if len(finite)
                else 0.5
            )
            mesh = axes[row_index, column].pcolormesh(
                x_edges,
                tau_edges,
                ratio.T,
                shading="auto",
                cmap="coolwarm",
                norm=TwoSlopeNorm(
                    vmin=1.0 - distance, vcenter=1.0, vmax=1.0 + distance
                ),
            )
            figure.colorbar(mesh, ax=axes[row_index, column], pad=0.01, label="ratio")
        for column in range(3):
            axes[row_index, column].set_xlabel(correlation_labels[name])
            axes[row_index, column].set_ylabel(
                r"$\log_{10}(\tau_{\rm f}/[\mathrm{fm}/c])$"
            )
    figure.suptitle(
        rf"All C/A declusterings: anti-$k_T$ R={radius:g}, {pair.pt_label(pt_low, pt_high)}"
        + f"\n{sample_label}; "
        + r"$\tau_{\rm f}=\hbar cE/Q^2$; "
        + "ratio cells require at least 10 no-prehydro entries",
        fontsize=11,
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    for suffix in ("pdf", "png"):
        figure.savefig(output_base.with_suffix(f".{suffix}"), dpi=180)
    plt.close(figure)


def run(args: argparse.Namespace) -> dict[str, Any]:
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    first_metadata = json.loads(args.metadata_alpha037.read_text())
    second_metadata = json.loads(args.metadata_alpha0335.read_text())
    for metadata, label in ((first_metadata, "alpha037"), (second_metadata, "alpha0335")):
        if int(metadata["pairCount"]) != args.expected_events:
            raise ValueError(f"{label} pair count is not {args.expected_events}")
    for field in ("pairCount", "weightSum", "sigmaMergedMb", "effectiveEventCount"):
        if not numeric_match(str(first_metadata[field]), str(second_metadata[field])):
            raise ValueError(f"pair metadata differs in {field}")

    spectrum_rows = merge_table(
        read_tsv(args.spectra_alpha037),
        read_tsv(args.spectra_alpha0335),
        SPECTRUM_STRUCTURE,
        (
            "weighted_yield_per_event",
            "weighted_yield_stat_error",
            "selected_jet_normalized_density",
            "selected_jet_normalized_density_stat_error",
            "differential_cross_section_mb",
            "cross_section_stat_error_mb",
            "raw_entries_in_bin",
            "raw_entries_total",
            "underflow_entries",
            "overflow_entries",
        ),
    )
    summary_rows = merge_table(
        read_tsv(args.summary_alpha037),
        read_tsv(args.summary_alpha0335),
        SUMMARY_STRUCTURE,
        tuple(
            field
            for field in read_tsv(args.summary_alpha037)[0]
            if field not in {"variant", *SUMMARY_STRUCTURE}
        ),
    )
    correlation_rows = merge_table(
        read_tsv(args.correlations_alpha037),
        read_tsv(args.correlations_alpha0335),
        CORRELATION_STRUCTURE,
        (
            "raw_entries",
            "selected_jet_normalized_density",
            "weighted_density_per_event",
            "weighted_density_mb",
        ),
    )
    for radius in pair.RADII:
        for pt_low, pt_high in pair.PT_INTERVALS:
            tag = f"R{int(round(10 * radius)):02d}_{pair.pt_tag(pt_low, pt_high)}"
            plot_spectra(
                spectrum_rows,
                radius,
                pt_low,
                pt_high,
                out_dir / f"{args.prefix}_{tag}_formation_time",
                args.sample_label,
            )
            plot_correlations(
                correlation_rows,
                radius,
                pt_low,
                pt_high,
                out_dir / f"{args.prefix}_{tag}_formation_time_correlations",
                args.sample_label,
            )
    write_tsv(out_dir / f"{args.prefix}_spectra.tsv", spectrum_rows)
    write_tsv(out_dir / f"{args.prefix}_summary.tsv", summary_rows)
    write_tsv(out_dir / f"{args.prefix}_correlations.tsv", correlation_rows)
    metadata = {
        "schemaVersion": "oo-v3-jet-formation-time-v3",
        "pairCount": args.expected_events,
        "strictNoPrehydroSpectrumAudit": "PASS",
        "strictNoPrehydroSummaryAudit": "PASS",
        "strictNoPrehydroCorrelationAudit": "PASS",
        "alpha037Metadata": str(args.metadata_alpha037.resolve()),
        "alpha0335Metadata": str(args.metadata_alpha0335.resolve()),
        "sampleLabel": args.sample_label,
        "formationTime": first_metadata["formationTime"],
        "caReclustering": (
            "FastJet cambridge_algorithm (p=0), E-scheme, Best strategy, "
            "R_CA=2*R_antiKt+1e-6"
        ),
        "singleSplitting": first_metadata["singleSplitting"],
        "negativeTreatment": first_metadata["negativeTreatment"],
        "interpretation": first_metadata["interpretation"],
    }
    (out_dir / f"{args.prefix}_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )
    return metadata


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    for kind in ("spectra", "summary", "correlations", "metadata"):
        result.add_argument(f"--{kind}-alpha037", required=True, type=Path)
        result.add_argument(f"--{kind}-alpha0335", required=True, type=Path)
    result.add_argument("--expected-events", required=True, type=int)
    result.add_argument("--out-dir", required=True, type=Path)
    result.add_argument("--prefix", default="oo5360_v3_jet_formation_time")
    result.add_argument(
        "--sample-label", default="same-seed V3 matched provisional sample"
    )
    return result


if __name__ == "__main__":
    raise SystemExit(0 if run(parser().parse_args()) else 1)
