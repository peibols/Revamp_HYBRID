#!/usr/bin/env python3
"""Plot event-weighted no/with-prehydro OO jet-variable comparisons."""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import math
from pathlib import Path

import awkward as ak
import numpy as np
import uproot


VARIANTS = ("noPrehydro", "withPrehydro")
RADIUS_DIGITS = (1, 2, 4, 8)
VARIANT_LABELS = {
    "noPrehydro": "No pre-hydro",
    "withPrehydro": "Plan B pre-hydro",
}
VARIANT_COLORS = {
    "noPrehydro": "#0072B2",
    "withPrehydro": "#D55E00",
}
VARIANT_MARKERS = {
    "noPrehydro": "o",
    "withPrehydro": "s",
}


@dataclasses.dataclass(frozen=True)
class VariableSpec:
    key: str
    branch_suffix: str
    label: str
    ylabel: str
    edges: np.ndarray
    group: str
    xscale: str = "linear"
    yscale: str = "linear"
    soft_drop_failure_bin: bool = False
    physical_ticks: tuple[float, ...] = ()


@dataclasses.dataclass
class HistogramResult:
    values: np.ndarray
    errors: np.ndarray
    weighted_matrix: np.ndarray
    raw_entries: int
    weighted_entries: float
    nonfinite_entries: int
    underflow_entries: int
    overflow_entries: int


@dataclasses.dataclass
class IntegratedResult:
    raw_jets: int
    events_with_jets: int
    weighted_jets: float
    cross_section: float
    error: float


@dataclasses.dataclass
class MeanResult:
    mean: float
    error: float
    raw_entries: int
    weighted_entries: float


def variable_specs(
    radius_digit: int, pt_min: float, pt_max: float | None = None
) -> list[VariableSpec]:
    if radius_digit == 1:
        mass_edges = np.array(
            [-32, -16, -8, -4, -2, -1, 0, 0.5, 1, 1.5, 2, 3, 4, 6, 10, 16, 32, 64],
            dtype=float,
        )
        multiplicity_edges = np.arange(0.5, 41.5, 1.0)
        rg_edges = np.linspace(0.0, 0.125, 21)
        girth_edges = np.linspace(0.0, 0.08, 21)
        max_kt_edges = np.geomspace(0.0005, 50.0, 20)
    elif radius_digit == 2:
        mass_edges = np.array(
            [-16, -8, -4, -2, 0, 1, 2, 3, 4, 6, 8, 12, 20, 32, 64, 128],
            dtype=float,
        )
        multiplicity_edges = np.arange(0.5, 61.5, 2.0)
        rg_edges = np.linspace(0.0, 0.25, 21)
        girth_edges = np.linspace(0.0, 0.16, 21)
        max_kt_edges = np.geomspace(0.0025, 100.0, 20)
    elif radius_digit == 4:
        mass_edges = np.array(
            [-32, -16, -8, -4, -2, 0, 2, 4, 6, 8, 10, 14, 20, 32, 64, 128, 256],
            dtype=float,
        )
        multiplicity_edges = np.arange(0.5, 91.5, 3.0)
        rg_edges = np.linspace(0.0, 0.5, 21)
        girth_edges = np.linspace(0.0, 0.30, 21)
        max_kt_edges = np.geomspace(0.005, 200.0, 20)
    elif radius_digit == 8:
        mass_edges = np.array(
            [-32, -16, -8, -4, -2, 0, 2, 4, 6, 8, 10, 14, 20, 28, 40, 64, 128, 256, 512],
            dtype=float,
        )
        multiplicity_edges = np.arange(0.5, 137.5, 4.0)
        rg_edges = np.linspace(0.0, 1.0, 21)
        girth_edges = np.linspace(0.0, 0.55, 23)
        max_kt_edges = np.geomspace(0.01, 400.0, 20)
    else:
        raise ValueError(f"unsupported radius digit {radius_digit}")

    rg_bin_width = float(rg_edges[1] - rg_edges[0])
    rg_edges = np.concatenate(([float(rg_edges[0] - rg_bin_width)], rg_edges))
    zg_edges = np.concatenate(([0.075], np.linspace(0.1, 0.5, 17)))

    if pt_max is None:
        pt_edges = np.geomspace(pt_min, 1600.0, 17)
        pt_xscale = "log"
    else:
        pt_edges = np.linspace(pt_min, pt_max, 11)
        pt_xscale = "linear"

    return [
        VariableSpec(
            "pt",
            "Pt",
            r"$p_T^{\mathrm{jet}}$ [GeV]",
            r"$d\sigma_{\mathrm{jet}}/dp_T$ [mb/GeV]",
            pt_edges,
            "kinematics",
            xscale=pt_xscale,
            yscale="log",
        ),
        VariableSpec(
            "eta",
            "Eta",
            r"$\eta_{\mathrm{jet}}$",
            r"$d\sigma_{\mathrm{jet}}/d\eta$ [mb]",
            np.linspace(-5.0, 5.0, 21),
            "kinematics",
        ),
        VariableSpec(
            "phi",
            "Phi",
            r"$\phi_{\mathrm{jet}}$",
            r"$d\sigma_{\mathrm{jet}}/d\phi$ [mb/rad]",
            np.linspace(-math.pi, math.pi, 19),
            "kinematics",
        ),
        VariableSpec(
            "mass",
            "mass",
            r"signed $m_{\mathrm{jet}}$ [GeV]",
            r"$d\sigma_{\mathrm{jet}}/dm$ [mb/GeV]",
            mass_edges,
            "kinematics",
            xscale="symlog",
            yscale="log",
        ),
        VariableSpec(
            "zg",
            "Zg",
            r"$z_g$",
            r"$d\sigma_{\mathrm{jet}}/dz_g$ [mb]",
            zg_edges,
            "substructure",
            soft_drop_failure_bin=True,
            physical_ticks=(0.2, 0.3, 0.4, 0.5),
        ),
        VariableSpec(
            "rg",
            "Rg",
            r"$R_g$",
            r"$d\sigma_{\mathrm{jet}}/dR_g$ [mb]",
            rg_edges,
            "substructure",
            soft_drop_failure_bin=True,
            physical_ticks=tuple(
                float(value) for value in np.linspace(0.0, rg_edges[-1], 6)[1:]
            ),
        ),
        VariableSpec(
            "mult",
            "Mult",
            r"constituent multiplicity",
            r"$d\sigma_{\mathrm{jet}}/dN$ [mb]",
            multiplicity_edges,
            "substructure",
            yscale="log",
        ),
        VariableSpec(
            "ptd",
            "PtD",
            r"$p_TD$",
            r"$d\sigma_{\mathrm{jet}}/d(p_TD)$ [mb]",
            np.linspace(0.0, 1.0, 21),
            "substructure",
        ),
        VariableSpec(
            "girth",
            "G",
            r"girth $g$",
            r"$d\sigma_{\mathrm{jet}}/dg$ [mb]",
            girth_edges,
            "substructure",
        ),
        VariableSpec(
            "maxkt",
            "MaxKt",
            r"maximum $k_T$ splitting [GeV]",
            r"$d\sigma_{\mathrm{jet}}/dk_T$ [mb/GeV]",
            max_kt_edges,
            "substructure",
            xscale="log",
            yscale="log",
        ),
    ]


def jet_pt_selection(
    pt_values: ak.Array | np.ndarray, pt_min: float, pt_max: float | None
) -> ak.Array | np.ndarray:
    selection = np.isfinite(pt_values) & (pt_values > pt_min)
    if pt_max is not None:
        selection = selection & (pt_values <= pt_max)
    return selection


def soft_drop_histogram_values(
    values: ak.Array, valid: ak.Array, spec: VariableSpec
) -> ak.Array:
    if not spec.soft_drop_failure_bin:
        raise ValueError(f"{spec.key} does not define a failed-Soft-Drop bin")
    failure_value = 0.5 * (spec.edges[0] + spec.edges[1])
    return ak.where(valid, values, failure_value)


def pt_range_label(pt_min: float, pt_max: float | None) -> str:
    if pt_max is None:
        return rf"$p_T^{{\rm jet}}>{pt_min:g}$ GeV"
    return rf"${pt_min:g}<p_T^{{\rm jet}}\leq {pt_max:g}$ GeV"


def jackknife_error(leave_one_out: np.ndarray) -> np.ndarray:
    if leave_one_out.ndim == 1:
        leave_one_out = leave_one_out[:, np.newaxis]
        squeeze = True
    else:
        squeeze = False
    finite = np.all(np.isfinite(leave_one_out), axis=0)
    errors = np.full(leave_one_out.shape[1], np.nan, dtype=float)
    if np.any(finite):
        values = leave_one_out[:, finite]
        means = np.mean(values, axis=0)
        differences = values - means
        count = values.shape[0]
        errors[finite] = np.sqrt((count - 1.0) / count * np.sum(differences**2, axis=0))
    return errors[0] if squeeze else errors


def pythia_parallel_factor(weights: np.ndarray, sigma_gen: np.ndarray) -> tuple[float, float, float]:
    weight_sum = float(np.sum(weights))
    weighted_sigma_sum = float(np.sum(weights * sigma_gen))
    if weight_sum <= 0.0 or weighted_sigma_sum <= 0.0:
        raise ValueError("invalid PYTHIA aggregate")
    sigma_merged = weighted_sigma_sum / weight_sum
    return weighted_sigma_sum / (weight_sum * weight_sum), weight_sum, sigma_merged


def histogram_matrix(
    values: np.ndarray,
    event_indices: np.ndarray,
    weights: np.ndarray,
    edges: np.ndarray,
) -> tuple[np.ndarray, int, int]:
    bin_count = len(edges) - 1
    indices = np.searchsorted(edges, values, side="right") - 1
    indices[values == edges[-1]] = bin_count - 1
    underflow = int(np.count_nonzero(indices < 0))
    overflow = int(np.count_nonzero(indices >= bin_count))
    accepted = (indices >= 0) & (indices < bin_count)
    raw_matrix = np.zeros((len(weights), bin_count), dtype=np.float64)
    np.add.at(raw_matrix, (event_indices[accepted], indices[accepted]), 1.0)
    return raw_matrix * weights[:, np.newaxis], underflow, overflow


def differential_histogram(
    jagged_values: ak.Array,
    weights: np.ndarray,
    sigma_gen: np.ndarray,
    edges: np.ndarray,
) -> HistogramResult:
    lengths_before_finite = ak.to_numpy(ak.num(jagged_values, axis=1))
    finite_mask = np.isfinite(jagged_values)
    finite_values = jagged_values[finite_mask]
    lengths = ak.to_numpy(ak.num(finite_values, axis=1))
    flat_values = ak.to_numpy(ak.flatten(finite_values))
    event_indices = np.repeat(np.arange(len(weights), dtype=np.int64), lengths)
    matrix, underflow, overflow = histogram_matrix(flat_values, event_indices, weights, edges)

    factor, weight_sum, _ = pythia_parallel_factor(weights, sigma_gen)
    widths = np.diff(edges)
    totals = np.sum(matrix, axis=0)
    values = factor * totals / widths

    weighted_sigma_sum = float(np.sum(weights * sigma_gen))
    leave_weight_sum = weight_sum - weights
    leave_weighted_sigma = weighted_sigma_sum - weights * sigma_gen
    leave_totals = totals[np.newaxis, :] - matrix
    leave_values = (
        leave_weighted_sigma[:, np.newaxis]
        * leave_totals
        / (leave_weight_sum[:, np.newaxis] ** 2)
        / widths[np.newaxis, :]
    )
    errors = jackknife_error(leave_values)
    return HistogramResult(
        values=values,
        errors=errors,
        weighted_matrix=matrix,
        raw_entries=int(np.sum(lengths)),
        weighted_entries=float(np.sum(weights * lengths)),
        nonfinite_entries=int(np.sum(lengths_before_finite - lengths)),
        underflow_entries=underflow,
        overflow_entries=overflow,
    )


def integrated_cross_section(
    selected_counts: np.ndarray,
    weights: np.ndarray,
    sigma_gen: np.ndarray,
) -> IntegratedResult:
    factor, weight_sum, _ = pythia_parallel_factor(weights, sigma_gen)
    weighted_counts = weights * selected_counts
    total = float(np.sum(weighted_counts))
    weighted_sigma_sum = float(np.sum(weights * sigma_gen))
    leave_values = (
        (weighted_sigma_sum - weights * sigma_gen)
        * (total - weighted_counts)
        / (weight_sum - weights) ** 2
    )
    return IntegratedResult(
        raw_jets=int(np.sum(selected_counts)),
        events_with_jets=int(np.count_nonzero(selected_counts)),
        weighted_jets=total,
        cross_section=factor * total,
        error=float(jackknife_error(leave_values)),
    )


def paired_ratio(
    numerator_matrix: np.ndarray,
    denominator_matrix: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    numerator = np.sum(numerator_matrix, axis=0)
    denominator = np.sum(denominator_matrix, axis=0)
    ratio = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan),
        where=denominator != 0.0,
    )
    leave_denominator = denominator[np.newaxis, :] - denominator_matrix
    leave_numerator = numerator[np.newaxis, :] - numerator_matrix
    leave_ratio = np.divide(
        leave_numerator,
        leave_denominator,
        out=np.full_like(leave_numerator, np.nan),
        where=leave_denominator != 0.0,
    )
    return ratio, jackknife_error(leave_ratio)


def paired_integrated_ratio(
    numerator_counts: np.ndarray,
    denominator_counts: np.ndarray,
    weights: np.ndarray,
) -> tuple[float, float]:
    numerator_by_event = weights * numerator_counts
    denominator_by_event = weights * denominator_counts
    numerator = float(np.sum(numerator_by_event))
    denominator = float(np.sum(denominator_by_event))
    ratio = numerator / denominator
    leave_ratio = (numerator - numerator_by_event) / (denominator - denominator_by_event)
    return ratio, float(jackknife_error(leave_ratio))


def weighted_mean(jagged_values: ak.Array, weights: np.ndarray) -> MeanResult:
    finite_values = jagged_values[np.isfinite(jagged_values)]
    lengths = ak.to_numpy(ak.num(finite_values, axis=1))
    sums = ak.to_numpy(ak.sum(finite_values, axis=1))
    weighted_counts = weights * lengths
    weighted_sums = weights * sums
    denominator = float(np.sum(weighted_counts))
    numerator = float(np.sum(weighted_sums))
    if denominator == 0.0:
        return MeanResult(math.nan, math.nan, 0, 0.0)
    leave_values = (numerator - weighted_sums) / (denominator - weighted_counts)
    return MeanResult(
        mean=numerator / denominator,
        error=float(jackknife_error(leave_values)),
        raw_entries=int(np.sum(lengths)),
        weighted_entries=denominator,
    )


def ratio_limits(ratio: np.ndarray, error: np.ndarray) -> tuple[float, float]:
    stable = np.isfinite(ratio) & np.isfinite(error) & (error < 0.35)
    if not np.any(stable):
        return 0.5, 1.5
    low = min(0.8, float(np.min(ratio[stable] - error[stable])))
    high = max(1.2, float(np.max(ratio[stable] + error[stable])))
    return max(0.25, low - 0.05), min(2.0, high + 0.05)


def configure_x_axis(axis, spec: VariableSpec) -> None:
    if spec.xscale == "log":
        axis.set_xscale("log")
    elif spec.xscale == "symlog":
        axis.set_xscale("symlog", linthresh=2.0, linscale=0.7)
    axis.set_xlim(spec.edges[0], spec.edges[-1])
    if spec.soft_drop_failure_bin:
        failure_center = 0.5 * (spec.edges[0] + spec.edges[1])
        axis.set_xticks([failure_center, *spec.physical_ticks])
        axis.set_xticklabels(
            ["SD fail", *(f"{value:g}" for value in spec.physical_ticks)]
        )


def draw_panel(
    figure,
    outer_spec,
    spec: VariableSpec,
    results: dict[str, HistogramResult],
    ratio: np.ndarray,
    ratio_error: np.ndarray,
    show_legend: bool,
) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpecFromSubplotSpec

    inner = GridSpecFromSubplotSpec(
        2,
        1,
        subplot_spec=outer_spec,
        height_ratios=[3.1, 1.0],
        hspace=0.05,
    )
    upper = figure.add_subplot(inner[0])
    lower = figure.add_subplot(inner[1], sharex=upper)
    centers = np.sqrt(spec.edges[:-1] * spec.edges[1:]) if spec.xscale == "log" else 0.5 * (
        spec.edges[:-1] + spec.edges[1:]
    )
    if spec.soft_drop_failure_bin:
        for axis in (upper, lower):
            axis.axvspan(spec.edges[0], spec.edges[1], color="0.90", zorder=0)
            axis.axvline(spec.edges[1], color="0.55", linewidth=0.8, linestyle=":")

    for variant in VARIANTS:
        result = results[variant]
        finite = np.isfinite(result.values) & np.isfinite(result.errors)
        if spec.yscale == "log":
            finite &= result.values > 0.0
        upper.stairs(
            result.values,
            spec.edges,
            color=VARIANT_COLORS[variant],
            linewidth=1.5,
            label=VARIANT_LABELS[variant],
        )
        upper.errorbar(
            centers[finite],
            result.values[finite],
            yerr=result.errors[finite],
            color=VARIANT_COLORS[variant],
            marker=VARIANT_MARKERS[variant],
            markersize=2.8,
            linestyle="none",
            linewidth=0.8,
            capsize=1.5,
        )
    if spec.yscale == "log":
        upper.set_yscale("log")
        positive = np.concatenate(
            [
                result.values[np.isfinite(result.values) & (result.values > 0.0)]
                for result in results.values()
            ]
        )
        if positive.size:
            upper.set_ylim(max(float(np.min(positive)) * 0.35, 1e-12), float(np.max(positive)) * 3.0)
    else:
        maximum = max(
            float(np.nanmax(result.values + np.nan_to_num(result.errors, nan=0.0)))
            for result in results.values()
        )
        upper.set_ylim(0.0, maximum * 1.28 if maximum > 0.0 else 1.0)
    upper.set_ylabel(spec.ylabel, fontsize=8.5)
    upper.tick_params(labelbottom=False, labelsize=8)
    upper.grid(alpha=0.22)
    if show_legend:
        upper.legend(frameon=False, fontsize=8, loc="best")

    finite_ratio = np.isfinite(ratio) & np.isfinite(ratio_error)
    lower.axhline(1.0, color="0.45", linewidth=0.9)
    lower.errorbar(
        centers[finite_ratio],
        ratio[finite_ratio],
        yerr=ratio_error[finite_ratio],
        color="#333333",
        marker="o",
        markersize=2.5,
        linestyle="none",
        linewidth=0.75,
        capsize=1.2,
    )
    lower.set_ylim(*ratio_limits(ratio, ratio_error))
    lower.set_ylabel("Pre/no", fontsize=8)
    lower.set_xlabel(spec.label, fontsize=9)
    lower.tick_params(labelsize=8)
    lower.grid(alpha=0.22)
    configure_x_axis(upper, spec)
    configure_x_axis(lower, spec)


def plot_group(
    out_dir: Path,
    prefix: str,
    radius_digit: int,
    pt_min: float,
    pt_max: float | None,
    specs: list[VariableSpec],
    histograms: dict[str, dict[str, HistogramResult]],
    ratios: dict[str, tuple[np.ndarray, np.ndarray]],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    group = specs[0].group
    if group == "kinematics":
        rows, columns, size = 2, 2, (10.8, 8.0)
    else:
        rows, columns, size = 2, 3, (15.2, 8.0)
    figure = plt.figure(figsize=size)
    grid = figure.add_gridspec(rows, columns, wspace=0.33, hspace=0.33)
    for index, spec in enumerate(specs):
        draw_panel(
            figure,
            grid[index // columns, index % columns],
            spec,
            histograms[spec.key],
            ratios[spec.key][0],
            ratios[spec.key][1],
            show_legend=index == 0,
        )
    radius = radius_digit / 10.0
    figure.suptitle(
        rf"O16+O16 5.36 TeV, 0--5% diagnostic; anti-$k_T$ $R={radius:.1f}$, "
        f"4MomSub {pt_range_label(pt_min, pt_max)}",
        fontsize=13,
        y=0.992,
    )
    figure.text(
        0.5,
        0.012,
        "PythiaParallel weighted cross sections; paired delete-one-run jackknife. "
        "No additional jet-eta cut. First Zg/Rg bin is SoftDropValid=0; physical bins are valid jets.",
        ha="center",
        fontsize=8.5,
    )
    figure.subplots_adjust(top=0.93, bottom=0.09)
    stem = out_dir / f"{prefix}_R0{radius_digit}_{group}"
    figure.savefig(stem.with_suffix(".pdf"))
    figure.savefig(stem.with_suffix(".png"), dpi=180)
    plt.close(figure)


def format_float(value: float) -> str:
    return "nan" if not math.isfinite(value) else f"{value:.12e}"


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--pt-min", type=float, default=30.0)
    parser.add_argument(
        "--pt-max",
        type=float,
        help="optional inclusive upper pT bound; the lower bound is exclusive",
    )
    parser.add_argument("--prefix", default="oo5360_jet_variables_pt30")
    return parser


def run(args: argparse.Namespace) -> dict[str, object]:
    if args.pt_min <= 0.0:
        raise ValueError("--pt-min must be positive")
    if args.pt_max is not None and args.pt_max <= args.pt_min:
        raise ValueError("--pt-max must be greater than --pt-min")
    input_root = args.input_root.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    all_suffixes = {
        spec.branch_suffix
        for radius in RADIUS_DIGITS
        for spec in variable_specs(radius, args.pt_min, args.pt_max)
    }
    all_suffixes.add("SoftDropValid")
    tree_arrays: dict[str, ak.Array] = {}
    with uproot.open(input_root) as root_file:
        pairs = root_file["Pairs"].arrays(["pairId", "eventWeight", "sigmaGen", "seed"], library="np")
        pair_ids = pairs["pairId"]
        weights = pairs["eventWeight"].astype(np.float64)
        sigma_gen = pairs["sigmaGen"].astype(np.float64)
        if len(pair_ids) == 0 or len(np.unique(pair_ids)) != len(pair_ids):
            raise ValueError("Pairs tree has empty or duplicate pairId values")
        if len(np.unique(pairs["seed"])) != len(pair_ids):
            raise ValueError("Pairs tree contains duplicate hard-event seeds")
        if not np.all(np.isfinite(weights)) or not np.all(weights > 0.0):
            raise ValueError("event weights must be finite and positive")
        if not np.all(np.isfinite(sigma_gen)) or not np.all(sigma_gen > 0.0):
            raise ValueError("sigmaGen values must be finite and positive")

        for variant in VARIANTS:
            branch_names = ["pairId", "eventWeight", "sigmaGen"]
            for radius_digit in RADIUS_DIGITS:
                branch_names.extend(f"jet{radius_digit}{suffix}" for suffix in sorted(all_suffixes))
            arrays = root_file[f"{variant}/Jets"].arrays(branch_names, library="ak")
            if not np.array_equal(ak.to_numpy(arrays.pairId), pair_ids):
                raise ValueError(f"{variant} pairId ordering differs from Pairs")
            if not np.array_equal(ak.to_numpy(arrays.eventWeight), weights):
                raise ValueError(f"{variant} eventWeight differs from Pairs")
            if not np.array_equal(ak.to_numpy(arrays.sigmaGen), sigma_gen):
                raise ValueError(f"{variant} sigmaGen differs from Pairs")
            tree_arrays[variant] = arrays

    factor, weight_sum, sigma_merged = pythia_parallel_factor(weights, sigma_gen)
    effective_events = float(weight_sum * weight_sum / np.sum(weights * weights))
    weighted_sigma_sum = float(np.sum(weights * sigma_gen))
    histogram_rows: list[list[object]] = []
    summary_rows: list[list[object]] = []
    moment_rows: list[list[object]] = []
    metadata: dict[str, object] = {
        "inputRoot": str(input_root),
        "inputBytes": input_root.stat().st_size,
        "pairCount": len(pair_ids),
        "uniqueSeedCount": len(np.unique(pairs["seed"])),
        "ptMinGeV": args.pt_min,
        "ptMaxGeV": args.pt_max,
        "ptSelection": (
            f"{args.pt_min:g} < corrected jet pT <= {args.pt_max:g} GeV"
            if args.pt_max is not None
            else f"corrected jet pT > {args.pt_min:g} GeV"
        ),
        "selectionJetPtBranch": "4MomSub-corrected jetPt",
        "additionalJetEtaCut": None,
        "weightSum": weight_sum,
        "weightedSigmaSum": weighted_sigma_sum,
        "sigmaMergedMb": sigma_merged,
        "crossSectionFactor": factor,
        "effectiveEventCountFromWeights": effective_events,
        "normalization": "dSigma/dx = sum(w*sigmaGen)/sum(w)^2 * sum(w*n_bin)/bin_width",
        "uncertainty": "paired delete-one-run jackknife",
        "softDropHistogramConvention": (
            "Zg and Rg bin 0 contains every selected SoftDropValid=0 jet at a finite "
            "sentinel; its width equals one physical bin, so height times width is the "
            "failed-jet cross section. Remaining bins contain SoftDropValid=1 jets."
        ),
        "softDropSelection": "all selected jets; bin 0 is SoftDropValid=0",
        "softDropMoments": "weighted Zg and Rg means use SoftDropValid=1 jets only",
        "radii": {},
    }

    for radius_digit in RADIUS_DIGITS:
        specs = variable_specs(radius_digit, args.pt_min, args.pt_max)
        selections: dict[str, ak.Array] = {}
        selected_counts: dict[str, np.ndarray] = {}
        integrated: dict[str, IntegratedResult] = {}
        for variant in VARIANTS:
            pt_values = tree_arrays[variant][f"jet{radius_digit}Pt"]
            selection = jet_pt_selection(pt_values, args.pt_min, args.pt_max)
            selections[variant] = selection
            counts = ak.to_numpy(ak.num(pt_values[selection], axis=1)).astype(np.float64)
            selected_counts[variant] = counts
            integrated[variant] = integrated_cross_section(counts, weights, sigma_gen)
        integrated_ratio, integrated_ratio_error = paired_integrated_ratio(
            selected_counts["withPrehydro"], selected_counts["noPrehydro"], weights
        )

        radius_metadata: dict[str, object] = {
            "integratedPreOverNoRatio": integrated_ratio,
            "integratedPreOverNoRatioError": integrated_ratio_error,
            "softDropFailureBins": {
                spec.key: {
                    "binIndex": 0,
                    "binLow": float(spec.edges[0]),
                    "binHigh": float(spec.edges[1]),
                    "binCenter": float(0.5 * (spec.edges[0] + spec.edges[1])),
                }
                for spec in specs
                if spec.soft_drop_failure_bin
            },
            "variants": {},
        }
        for variant in VARIANTS:
            soft_drop_valid = (
                tree_arrays[variant][f"jet{radius_digit}SoftDropValid"][selections[variant]]
                == 1
            )
            valid_counts = ak.to_numpy(ak.sum(soft_drop_valid, axis=1)).astype(np.float64)
            selected = integrated[variant]
            valid_weighted = float(np.sum(weights * valid_counts))
            soft_drop_fraction = valid_weighted / selected.weighted_jets
            failed_counts = selected_counts[variant] - valid_counts
            failed_weighted = float(np.sum(weights * failed_counts))
            failed_fraction = failed_weighted / selected.weighted_jets
            radius_metadata["variants"][variant] = {
                "rawSelectedJets": selected.raw_jets,
                "eventsWithSelectedJets": selected.events_with_jets,
                "weightedSelectedJets": selected.weighted_jets,
                "integratedJetCrossSectionMb": selected.cross_section,
                "integratedJetCrossSectionErrorMb": selected.error,
                "rawSoftDropValidJets": int(np.sum(valid_counts)),
                "weightedSoftDropValidFraction": soft_drop_fraction,
                "rawSoftDropFailedJets": int(np.sum(failed_counts)),
                "weightedSoftDropFailedFraction": failed_fraction,
            }
            summary_rows.append(
                [
                    f"0.{radius_digit}",
                    variant,
                    selected.raw_jets,
                    selected.events_with_jets,
                    format_float(selected.weighted_jets),
                    format_float(selected.cross_section),
                    format_float(selected.error),
                    format_float(integrated_ratio),
                    format_float(integrated_ratio_error),
                    int(np.sum(valid_counts)),
                    format_float(soft_drop_fraction),
                    int(np.sum(failed_counts)),
                    format_float(failed_fraction),
                    len(pair_ids),
                    format_float(weight_sum),
                    format_float(sigma_merged),
                    format_float(effective_events),
                ]
            )
        metadata["radii"][f"0.{radius_digit}"] = radius_metadata

        histograms: dict[str, dict[str, HistogramResult]] = {}
        ratios: dict[str, tuple[np.ndarray, np.ndarray]] = {}
        for spec in specs:
            histograms[spec.key] = {}
            jagged_by_variant: dict[str, ak.Array] = {}
            for variant in VARIANTS:
                values = tree_arrays[variant][f"jet{radius_digit}{spec.branch_suffix}"][
                    selections[variant]
                ]
                histogram_values = values
                moment_values = values
                if spec.soft_drop_failure_bin:
                    valid = (
                        tree_arrays[variant][f"jet{radius_digit}SoftDropValid"][
                            selections[variant]
                        ]
                        == 1
                    )
                    histogram_values = soft_drop_histogram_values(values, valid, spec)
                    moment_values = values[valid]
                jagged_by_variant[variant] = moment_values
                result = differential_histogram(
                    histogram_values, weights, sigma_gen, spec.edges
                )
                if spec.soft_drop_failure_bin:
                    failed_counts = ak.to_numpy(ak.sum(~valid, axis=1)).astype(np.float64)
                    expected_failure_bin = weights * failed_counts
                    if (
                        result.raw_entries != integrated[variant].raw_jets
                        or result.nonfinite_entries != 0
                        or result.underflow_entries != 0
                        or result.overflow_entries != 0
                        or not np.array_equal(
                            result.weighted_matrix[:, 0], expected_failure_bin
                        )
                    ):
                        raise ValueError(
                            f"R=0.{radius_digit} {variant} {spec.key} failed-bin "
                            "closure does not reproduce every selected jet"
                        )
                histograms[spec.key][variant] = result
            ratios[spec.key] = paired_ratio(
                histograms[spec.key]["withPrehydro"].weighted_matrix,
                histograms[spec.key]["noPrehydro"].weighted_matrix,
            )

            centers = (
                np.sqrt(spec.edges[:-1] * spec.edges[1:])
                if spec.xscale == "log"
                else 0.5 * (spec.edges[:-1] + spec.edges[1:])
            )
            ratio, ratio_error = ratios[spec.key]
            for variant in VARIANTS:
                result = histograms[spec.key][variant]
                for bin_index, (low, high, center, value, error) in enumerate(
                    zip(spec.edges[:-1], spec.edges[1:], centers, result.values, result.errors)
                ):
                    histogram_rows.append(
                        [
                            f"0.{radius_digit}",
                            spec.group,
                            spec.key,
                            variant,
                            bin_index,
                            (
                                "softdrop_failed"
                                if spec.soft_drop_failure_bin and bin_index == 0
                                else "physical"
                            ),
                            format_float(float(low)),
                            format_float(float(high)),
                            format_float(float(center)),
                            format_float(float(value)),
                            format_float(float(error)),
                            format_float(float(ratio[bin_index])),
                            format_float(float(ratio_error[bin_index])),
                            result.raw_entries,
                            format_float(result.weighted_entries),
                            result.nonfinite_entries,
                            result.underflow_entries,
                            result.overflow_entries,
                        ]
                    )

            means = {variant: weighted_mean(jagged_by_variant[variant], weights) for variant in VARIANTS}
            # A paired mean-ratio jackknife is reconstructed from event-level sums and counts.
            event_sums: dict[str, np.ndarray] = {}
            event_counts: dict[str, np.ndarray] = {}
            for variant in VARIANTS:
                finite = jagged_by_variant[variant][np.isfinite(jagged_by_variant[variant])]
                event_sums[variant] = weights * ak.to_numpy(ak.sum(finite, axis=1))
                event_counts[variant] = weights * ak.to_numpy(ak.num(finite, axis=1))
            totals = {variant: float(np.sum(event_sums[variant])) for variant in VARIANTS}
            counts = {variant: float(np.sum(event_counts[variant])) for variant in VARIANTS}
            leave_no = (totals["noPrehydro"] - event_sums["noPrehydro"]) / (
                counts["noPrehydro"] - event_counts["noPrehydro"]
            )
            leave_pre = (totals["withPrehydro"] - event_sums["withPrehydro"]) / (
                counts["withPrehydro"] - event_counts["withPrehydro"]
            )
            mean_difference = means["withPrehydro"].mean - means["noPrehydro"].mean
            mean_difference_error = float(jackknife_error(leave_pre - leave_no))
            if spec.key in {"eta", "phi"} or means["noPrehydro"].mean == 0.0:
                mean_ratio = math.nan
                mean_ratio_error = math.nan
            else:
                mean_ratio = means["withPrehydro"].mean / means["noPrehydro"].mean
                mean_ratio_error = float(jackknife_error(leave_pre / leave_no))
            for variant in VARIANTS:
                mean = means[variant]
                moment_rows.append(
                    [
                        f"0.{radius_digit}",
                        spec.key,
                        variant,
                        format_float(mean.mean),
                        format_float(mean.error),
                        mean.raw_entries,
                        format_float(mean.weighted_entries),
                        format_float(mean_ratio),
                        format_float(mean_ratio_error),
                        format_float(mean_difference),
                        format_float(mean_difference_error),
                    ]
                )

        for group in ("kinematics", "substructure"):
            plot_group(
                out_dir,
                args.prefix,
                radius_digit,
                args.pt_min,
                args.pt_max,
                [spec for spec in specs if spec.group == group],
                histograms,
                ratios,
            )

    with (out_dir / f"{args.prefix}_histograms.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "radius",
                "group",
                "variable",
                "variant",
                "bin_index",
                "bin_kind",
                "bin_low",
                "bin_high",
                "bin_center",
                "differential_cross_section_mb",
                "stat_error_mb",
                "pre_over_no_ratio",
                "ratio_stat_error",
                "raw_variable_entries",
                "weighted_variable_entries",
                "nonfinite_entries",
                "underflow_entries",
                "overflow_entries",
            ]
        )
        writer.writerows(histogram_rows)

    with (out_dir / f"{args.prefix}_summary.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "radius",
                "variant",
                "raw_selected_jets",
                "events_with_selected_jets",
                "weighted_selected_jets",
                "integrated_jet_cross_section_mb",
                "integrated_jet_cross_section_stat_error_mb",
                "pre_over_no_integrated_ratio",
                "ratio_stat_error",
                "raw_softdrop_valid_jets",
                "weighted_softdrop_valid_fraction",
                "raw_softdrop_failed_jets",
                "weighted_softdrop_failed_fraction",
                "pair_count",
                "weight_sum",
                "sigma_merged_mb",
                "effective_event_count",
            ]
        )
        writer.writerows(summary_rows)

    with (out_dir / f"{args.prefix}_moments.tsv").open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "radius",
                "variable",
                "variant",
                "weighted_mean",
                "mean_stat_error",
                "raw_entries",
                "weighted_entries",
                "pre_over_no_mean_ratio",
                "mean_ratio_stat_error",
                "pre_minus_no_mean",
                "mean_difference_stat_error",
            ]
        )
        writer.writerows(moment_rows)

    metadata_path = out_dir / f"{args.prefix}_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return metadata


def main() -> int:
    args = make_parser().parse_args()
    metadata = run(args)
    print(json.dumps(metadata, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
