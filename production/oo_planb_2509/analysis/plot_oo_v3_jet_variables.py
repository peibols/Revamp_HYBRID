#!/usr/bin/env python3
"""Plot strict three-way V3 OO jet-variable spectra from aligned pair tables."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
import numpy as np

import convert_oo_paired_to_root as paired
import plot_oo_jet_variables as base


VARIANTS = ("no_prehydro", "prehydro_alpha037", "prehydro_alpha0335")
RATIO_VARIANTS = ("prehydro_alpha037", "prehydro_alpha0335")
STYLES = {
    "no_prehydro": {
        "label": "No pre-hydro",
        "color": "#111111",
        "marker": None,
    },
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
STRUCTURAL_FIELDS = (
    "radius",
    "group",
    "variable",
    "bin_index",
    "bin_kind",
    "bin_low",
    "bin_high",
    "bin_center",
)
NO_SPECTRUM_FIELDS = (
    *STRUCTURAL_FIELDS,
    "differential_cross_section_mb",
    "stat_error_mb",
    "normalized_density",
    "normalized_density_stat_error",
    "raw_variable_entries",
    "weighted_variable_entries",
    "nonfinite_entries",
    "underflow_entries",
    "overflow_entries",
)
SUMMARY_NO_FIELDS = (
    "radius",
    "raw_selected_jets",
    "events_with_selected_jets",
    "weighted_selected_jets",
    "integrated_jet_cross_section_mb",
    "integrated_jet_cross_section_stat_error_mb",
    "raw_softdrop_valid_jets",
    "weighted_softdrop_valid_fraction",
    "raw_softdrop_failed_jets",
    "weighted_softdrop_failed_fraction",
    "pair_count",
    "weight_sum",
    "sigma_merged_mb",
    "effective_event_count",
)
SUMMARY_COMMON_FIELDS = (
    "radius",
    "pair_count",
    "weight_sum",
    "sigma_merged_mb",
    "effective_event_count",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"empty TSV: {path}")
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty TSV: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def histogram_key(row: dict[str, str]) -> tuple[str, str, str, int]:
    return (
        row["radius"],
        row["group"],
        row["variable"],
        int(row["bin_index"]),
    )


def variant_map(
    rows: list[dict[str, str]], variant: str
) -> dict[tuple[str, str, str, int], dict[str, str]]:
    selected = {histogram_key(row): row for row in rows if row["variant"] == variant}
    if not selected:
        raise ValueError(f"no histogram rows for variant {variant}")
    return selected


def values_match(field: str, first: str, second: str) -> bool:
    if field in {
        "radius",
        "bin_low",
        "bin_high",
        "bin_center",
        "differential_cross_section_mb",
        "stat_error_mb",
        "normalized_density",
        "normalized_density_stat_error",
        "weighted_variable_entries",
        "integrated_jet_cross_section_mb",
        "integrated_jet_cross_section_stat_error_mb",
        "weighted_selected_jets",
        "weighted_softdrop_valid_fraction",
        "weighted_softdrop_failed_fraction",
        "weight_sum",
        "sigma_merged_mb",
        "effective_event_count",
    }:
        return math.isclose(
            float(first), float(second), rel_tol=1.0e-12, abs_tol=1.0e-15
        )
    return first == second


def assert_fields_match(
    first: dict[str, str],
    second: dict[str, str],
    fields: tuple[str, ...],
    description: str,
) -> None:
    for field in fields:
        if not values_match(field, first[field], second[field]):
            raise ValueError(
                f"{description} differs in {field}: {first[field]} != {second[field]}"
            )


def merge_histograms(
    alpha037_path: Path, alpha0335_path: Path
) -> list[dict[str, str]]:
    first = read_tsv(alpha037_path)
    second = read_tsv(alpha0335_path)
    maps = {
        "no037": variant_map(first, "noPrehydro"),
        "pre037": variant_map(first, "withPrehydro"),
        "no335": variant_map(second, "noPrehydro"),
        "pre335": variant_map(second, "withPrehydro"),
    }
    keys = set(maps["no037"])
    if any(set(rows) != keys for rows in maps.values()):
        raise ValueError("V3 histogram tables have different bin keys")

    output: list[dict[str, str]] = []
    for key in sorted(keys, key=lambda item: (float(item[0]), item[1], item[2], item[3])):
        no037 = maps["no037"][key]
        no335 = maps["no335"][key]
        assert_fields_match(no037, no335, NO_SPECTRUM_FIELDS, f"no-prehydro {key}")
        for candidate in (maps["pre037"][key], maps["pre335"][key]):
            assert_fields_match(no037, candidate, STRUCTURAL_FIELDS, f"bin {key}")
        for variant, source in (
            ("no_prehydro", no037),
            ("prehydro_alpha037", maps["pre037"][key]),
            ("prehydro_alpha0335", maps["pre335"][key]),
        ):
            converted = dict(source)
            converted["variant"] = variant
            if variant == "no_prehydro":
                converted["pre_over_no_ratio"] = "1.000000000000e+00"
                converted["ratio_stat_error"] = "0.000000000000e+00"
            output.append(converted)
    return output


def summary_map(
    rows: list[dict[str, str]], variant: str
) -> dict[str, dict[str, str]]:
    selected = {row["radius"]: row for row in rows if row["variant"] == variant}
    if set(selected) != {"0.1", "0.2", "0.4", "0.8"}:
        raise ValueError(f"incomplete summary radii for {variant}")
    return selected


def merge_summaries(
    alpha037_path: Path, alpha0335_path: Path, expected_events: int
) -> list[dict[str, str]]:
    first = read_tsv(alpha037_path)
    second = read_tsv(alpha0335_path)
    maps = {
        "no037": summary_map(first, "noPrehydro"),
        "pre037": summary_map(first, "withPrehydro"),
        "no335": summary_map(second, "noPrehydro"),
        "pre335": summary_map(second, "withPrehydro"),
    }
    output: list[dict[str, str]] = []
    for radius in ("0.1", "0.2", "0.4", "0.8"):
        no037 = maps["no037"][radius]
        no335 = maps["no335"][radius]
        assert_fields_match(
            no037, no335, SUMMARY_NO_FIELDS, f"no-prehydro summary R={radius}"
        )
        if int(no037["pair_count"]) != expected_events:
            raise ValueError(
                f"R={radius} pair_count={no037['pair_count']}, expected {expected_events}"
            )
        for candidate in (maps["pre037"][radius], maps["pre335"][radius]):
            assert_fields_match(
                no037,
                candidate,
                SUMMARY_COMMON_FIELDS,
                f"common summary metadata R={radius}",
            )
        for variant, source in (
            ("no_prehydro", no037),
            ("prehydro_alpha037", maps["pre037"][radius]),
            ("prehydro_alpha0335", maps["pre335"][radius]),
        ):
            if int(source["pair_count"]) != expected_events:
                raise ValueError(f"R={radius} {variant} has the wrong pair count")
            converted = dict(source)
            converted["variant"] = variant
            if variant == "no_prehydro":
                converted["pre_over_no_integrated_ratio"] = "1.000000000000e+00"
                converted["ratio_stat_error"] = "0.000000000000e+00"
            output.append(converted)
    return output


def rows_for_spec(
    rows: list[dict[str, str]], radius: float, spec: base.VariableSpec, variant: str
) -> list[dict[str, str]]:
    selected = [
        row
        for row in rows
        if row["variant"] == variant
        and row["variable"] == spec.key
        and math.isclose(float(row["radius"]), radius)
    ]
    selected.sort(key=lambda row: int(row["bin_index"]))
    if len(selected) != len(spec.edges) - 1:
        raise ValueError(
            f"R={radius} {spec.key} {variant} has {len(selected)} bins, "
            f"expected {len(spec.edges) - 1}"
        )
    edges = np.array(
        [float(selected[0]["bin_low"]), *[float(row["bin_high"]) for row in selected]]
    )
    if not np.allclose(edges, spec.edges, rtol=1.0e-12, atol=1.0e-14):
        raise ValueError(f"R={radius} {spec.key} bin edges disagree with the plot schema")
    return selected


def draw_panel(
    figure,
    outer_spec,
    spec: base.VariableSpec,
    rows: list[dict[str, str]],
    radius: float,
    show_legend: bool,
) -> None:
    inner = GridSpecFromSubplotSpec(
        2,
        1,
        subplot_spec=outer_spec,
        height_ratios=[3.1, 1.0],
        hspace=0.05,
    )
    upper = figure.add_subplot(inner[0])
    lower = figure.add_subplot(inner[1], sharex=upper)
    centers = (
        np.sqrt(spec.edges[:-1] * spec.edges[1:])
        if spec.xscale == "log"
        else 0.5 * (spec.edges[:-1] + spec.edges[1:])
    )
    if spec.soft_drop_failure_bin:
        for axis in (upper, lower):
            axis.axvspan(spec.edges[0], spec.edges[1], color="0.90", zorder=0)
            axis.axvline(spec.edges[1], color="0.55", linewidth=0.8, linestyle=":")

    selected: dict[str, list[dict[str, str]]] = {}
    positive: list[np.ndarray] = []
    for variant in VARIANTS:
        selected[variant] = rows_for_spec(rows, radius, spec, variant)
        values = np.array(
            [float(row["normalized_density"]) for row in selected[variant]]
        )
        errors = np.array(
            [float(row["normalized_density_stat_error"]) for row in selected[variant]]
        )
        finite = np.isfinite(values) & np.isfinite(errors)
        if spec.yscale == "log":
            finite &= values > 0.0
            positive.append(values[np.isfinite(values) & (values > 0.0)])
        style = STYLES[variant]
        if variant == "no_prehydro":
            upper.stairs(
                values,
                spec.edges,
                color=style["color"],
                linewidth=1.6,
                label=style["label"],
            )
        else:
            upper.errorbar(
                centers[finite],
                values[finite],
                yerr=errors[finite],
                color=style["color"],
                marker=style["marker"],
                markersize=2.8,
                linestyle="none",
                linewidth=0.8,
                capsize=1.3,
                label=style["label"],
            )

    if spec.yscale == "log":
        upper.set_yscale("log")
        finite_positive = np.concatenate(positive) if positive else np.array([])
        if finite_positive.size:
            upper.set_ylim(
                max(float(np.min(finite_positive)) * 0.35, 1.0e-12),
                float(np.max(finite_positive)) * 3.0,
            )
    else:
        maximum = 0.0
        for variant in VARIANTS:
            values = np.array(
                [
                    float(row["normalized_density"])
                    for row in selected[variant]
                ]
            )
            errors = np.array(
                [float(row["normalized_density_stat_error"]) for row in selected[variant]]
            )
            maximum = max(maximum, float(np.nanmax(values + np.nan_to_num(errors))))
        upper.set_ylim(0.0, maximum * 1.28 if maximum > 0.0 else 1.0)
    upper.set_ylabel(r"$(1/\sigma_{\rm jet})\,d\sigma/dx$", fontsize=8.2)
    upper.tick_params(labelbottom=False, labelsize=7.7)
    upper.grid(alpha=0.22)
    if show_legend:
        upper.legend(frameon=False, fontsize=7.5, loc="best")

    ratio_values: list[np.ndarray] = []
    ratio_errors: list[np.ndarray] = []
    lower.axhline(1.0, color="0.45", linewidth=0.9)
    for variant in RATIO_VARIANTS:
        ratio = np.array(
            [float(row["pre_over_no_ratio"]) for row in selected[variant]]
        )
        error = np.array([float(row["ratio_stat_error"]) for row in selected[variant]])
        ratio_values.append(ratio)
        ratio_errors.append(error)
        finite = np.isfinite(ratio) & np.isfinite(error)
        style = STYLES[variant]
        lower.errorbar(
            centers[finite],
            ratio[finite],
            yerr=error[finite],
            color=style["color"],
            marker=style["marker"],
            markersize=2.4,
            linestyle="none",
            linewidth=0.7,
            capsize=1.1,
        )
    lower.set_ylim(
        *base.ratio_limits(np.concatenate(ratio_values), np.concatenate(ratio_errors))
    )
    lower.set_ylabel("Pre/no", fontsize=7.6)
    lower.set_xlabel(spec.label, fontsize=8.7)
    lower.tick_params(labelsize=7.6)
    lower.grid(alpha=0.22)
    base.configure_x_axis(upper, spec)
    base.configure_x_axis(lower, spec)


def plot_groups(
    *,
    rows: list[dict[str, str]],
    out_dir: Path,
    prefix: str,
    pt_min: float,
    pt_max: float | None,
) -> list[Path]:
    plots: list[Path] = []
    for radius_digit in base.RADIUS_DIGITS:
        radius = radius_digit / 10.0
        specs = base.variable_specs(radius_digit, pt_min, pt_max)
        for group in ("kinematics", "substructure"):
            grouped = [spec for spec in specs if spec.group == group]
            if group == "kinematics":
                rows_count, columns, size = 2, 2, (10.8, 8.0)
            else:
                rows_count, columns, size = 2, 4, (18.0, 8.0)
            figure = plt.figure(figsize=size)
            grid = figure.add_gridspec(
                rows_count, columns, wspace=0.33, hspace=0.33
            )
            for index, spec in enumerate(grouped):
                draw_panel(
                    figure,
                    grid[index // columns, index % columns],
                    spec,
                    rows,
                    radius,
                    show_legend=index == 0,
                )
            figure.suptitle(
                rf"O+O 5.36 TeV, 0--5%, V3 same-seed matched; anti-$k_T$ "
                rf"$R={radius:.1f}$, 4MomSub {base.pt_range_label(pt_min, pt_max)}",
                fontsize=13,
                y=0.992,
            )
            figure.text(
                0.5,
                0.012,
                r"Per-variant $(1/\sigma_{\rm jet})d\sigma/dx$; paired delete-one-event jackknife. "
                "No additional jet-eta cut. First Zg/Rg bin is SoftDropValid=0.",
                ha="center",
                fontsize=8.3,
            )
            figure.subplots_adjust(top=0.93, bottom=0.09)
            path = out_dir / f"{prefix}_R0{radius_digit}_{group}.pdf"
            figure.savefig(path)
            figure.savefig(path.with_suffix(".png"), dpi=180)
            plt.close(figure)
            plots.append(path)
    return plots


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hist-alpha037", type=Path, required=True)
    parser.add_argument("--hist-alpha0335", type=Path, required=True)
    parser.add_argument("--summary-alpha037", type=Path, required=True)
    parser.add_argument("--summary-alpha0335", type=Path, required=True)
    parser.add_argument("--pt-min", type=float, required=True)
    parser.add_argument("--pt-max", type=float)
    parser.add_argument("--expected-events", type=int, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    if args.pt_min <= 0.0:
        raise ValueError("--pt-min must be positive")
    if args.pt_max is not None and args.pt_max <= args.pt_min:
        raise ValueError("--pt-max must exceed --pt-min")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    histograms = merge_histograms(args.hist_alpha037, args.hist_alpha0335)
    summaries = merge_summaries(
        args.summary_alpha037, args.summary_alpha0335, args.expected_events
    )
    histogram_output = args.out_dir / f"{args.prefix}_histograms.tsv"
    summary_output = args.out_dir / f"{args.prefix}_summary.tsv"
    write_tsv(histogram_output, histograms)
    write_tsv(summary_output, summaries)
    plots = plot_groups(
        rows=histograms,
        out_dir=args.out_dir,
        prefix=args.prefix,
        pt_min=args.pt_min,
        pt_max=args.pt_max,
    )

    radii: dict[str, object] = {}
    for radius in ("0.1", "0.2", "0.4", "0.8"):
        selected = {row["variant"]: row for row in summaries if row["radius"] == radius}
        radii[radius] = {
            "variants": {
                variant: {
                    "integratedJetCrossSectionMb": float(
                        selected[variant]["integrated_jet_cross_section_mb"]
                    ),
                    "integratedJetCrossSectionErrorMb": float(
                        selected[variant]["integrated_jet_cross_section_stat_error_mb"]
                    ),
                    "ratioToNoPrehydro": float(
                        selected[variant]["pre_over_no_integrated_ratio"]
                    ),
                    "ratioToNoPrehydroError": float(
                        selected[variant]["ratio_stat_error"]
                    ),
                }
                for variant in VARIANTS
            }
        }
    metadata = {
        "status": "PASS",
        "publicationStatus": "PROVISIONAL_COMPLETION_ORDER_MATCHED_DIAGNOSTIC",
        "expectedEvents": args.expected_events,
        "ptMinGeV": args.pt_min,
        "ptMaxGeV": args.pt_max,
        "selection": base.pt_range_label(args.pt_min, args.pt_max),
        "normalization": "PythiaParallel sigmaGen/sum(weight), applied once",
        "uncertainty": "paired delete-one-AA-event jackknife",
        "inputs": {
            "histAlpha037": str(args.hist_alpha037.resolve()),
            "histAlpha037Sha256": paired.sha256(args.hist_alpha037),
            "histAlpha0335": str(args.hist_alpha0335.resolve()),
            "histAlpha0335Sha256": paired.sha256(args.hist_alpha0335),
            "summaryAlpha037": str(args.summary_alpha037.resolve()),
            "summaryAlpha037Sha256": paired.sha256(args.summary_alpha037),
            "summaryAlpha0335": str(args.summary_alpha0335.resolve()),
            "summaryAlpha0335Sha256": paired.sha256(args.summary_alpha0335),
        },
        "radii": radii,
        "outputs": {
            "histograms": str(histogram_output),
            "summary": str(summary_output),
            "plots": [str(path) for path in plots],
        },
    }
    metadata_path = args.out_dir / f"{args.prefix}_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metadata, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
