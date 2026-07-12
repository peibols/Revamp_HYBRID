#!/usr/bin/env python3
"""Measure paired pre-hydro jet response versus effective-charge proxies."""

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


RADIUS_DIGITS = (2, 4, 8)
FLAVORS = ("all", "quark", "gluon")
FLAVOR_LABELS = {
    "all": "All jets",
    "quark": "Quark tagged",
    "gluon": "Gluon tagged",
}
FLAVOR_COLORS = {
    "all": "#222222",
    "quark": "#0072B2",
    "gluon": "#D55E00",
}
FLAVOR_MARKERS = {"all": "o", "quark": "s", "gluon": "^"}
CATEGORY_LABELS = (
    "single-like\n" + r"$N_{\rm eff}<3$",
    "intermediate\n" + r"$3\leq N_{\rm eff}<8$",
    "many-like\n" + r"$N_{\rm eff}\geq8$",
)
RADIUS_COLORS = {2: "#009E73", 4: "#0072B2", 8: "#D55E00"}
RADIUS_MARKERS = {2: "o", 4: "s", 8: "^"}


@dataclasses.dataclass(frozen=True)
class ProxySpec:
    key: str
    label: str
    edges: np.ndarray
    bin_labels: tuple[str, ...]


@dataclasses.dataclass
class BinnedResponse:
    mean: np.ndarray
    mean_error: np.ndarray
    pt_weighted_fraction: np.ndarray
    pt_weighted_fraction_error: np.ndarray
    raw_entries: np.ndarray
    weighted_entries: np.ndarray
    leave_mean: np.ndarray
    leave_pt_weighted_fraction: np.ndarray
    underflow_entries: int
    overflow_entries: int


def proxy_specs() -> tuple[ProxySpec, ...]:
    return (
        ProxySpec(
            "normal_neff",
            r"wake-excluded $N_{\rm eff}=1/(p_{TD}^{\rm normal})^2$",
            np.array([0.5, 2.0, 3.0, 5.0, 8.0, 13.0, 21.0, 34.0, 55.0]),
            ("<2", "2-3", "3-5", "5-8", "8-13", "13-21", "21-34", "34-55"),
        ),
        ProxySpec(
            "leading_normal_fraction",
            r"leading normal-hadron fraction $z_{\rm lead}^{\rm normal}$",
            np.linspace(0.0, 1.0, 11),
            tuple(f"{low:.1f}-{high:.1f}" for low, high in zip(np.linspace(0, 0.9, 10), np.linspace(0.1, 1, 10))),
        ),
        ProxySpec(
            "resolved_primary_prongs",
            r"hardest-branch count $n_{\rm SD}+1$",
            np.arange(0.5, 9.5, 1.0),
            tuple(str(value) for value in range(1, 9)),
        ),
        ProxySpec(
            "max_kt",
            r"largest C/A-tree $k_T$ splitting [GeV]",
            np.array([0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0]),
            ("<1", "1-2", "2-4", "4-8", "8-16", "16-32", "32-64", "64-128", "128-256", "256-512"),
        ),
    )


def jackknife_error(leave_one_out: np.ndarray) -> np.ndarray | float:
    values = np.asarray(leave_one_out, dtype=float)
    squeeze = values.ndim == 1
    if squeeze:
        values = values[:, np.newaxis]
    errors = np.full(values.shape[1], np.nan, dtype=float)
    for column in range(values.shape[1]):
        sample = values[:, column]
        if not np.all(np.isfinite(sample)) or len(sample) < 2:
            continue
        mean = float(np.mean(sample))
        errors[column] = math.sqrt(
            (len(sample) - 1.0) / len(sample) * float(np.sum((sample - mean) ** 2))
        )
    return float(errors[0]) if squeeze else errors


def flavor_mask(
    hard_parton_ids: np.ndarray,
    flavor: str,
    other_hard_parton_ids: np.ndarray | None = None,
) -> np.ndarray:
    absolute = np.abs(hard_parton_ids)
    if flavor == "all":
        return np.ones(len(hard_parton_ids), dtype=bool)
    if flavor == "quark":
        selected = (absolute >= 1) & (absolute <= 6)
    elif flavor == "gluon":
        selected = hard_parton_ids == 21
    else:
        raise ValueError(f"unsupported flavor {flavor}")
    if other_hard_parton_ids is not None:
        selected &= other_hard_parton_ids == hard_parton_ids
    return selected


def binned_response(
    proxy: np.ndarray,
    fractional_shift: np.ndarray,
    baseline_pt: np.ndarray,
    event_indices: np.ndarray,
    event_weights: np.ndarray,
    edges: np.ndarray,
) -> BinnedResponse:
    proxy = np.asarray(proxy, dtype=float)
    fractional_shift = np.asarray(fractional_shift, dtype=float)
    baseline_pt = np.asarray(baseline_pt, dtype=float)
    event_indices = np.asarray(event_indices, dtype=np.int64)
    if not (
        len(proxy) == len(fractional_shift) == len(baseline_pt) == len(event_indices)
    ):
        raise ValueError("jet arrays have inconsistent lengths")
    bin_count = len(edges) - 1
    finite = (
        np.isfinite(proxy)
        & np.isfinite(fractional_shift)
        & np.isfinite(baseline_pt)
        & (baseline_pt > 0.0)
    )
    indices = np.searchsorted(edges, proxy, side="right") - 1
    indices[proxy == edges[-1]] = bin_count - 1
    underflow = int(np.count_nonzero(finite & (indices < 0)))
    overflow = int(np.count_nonzero(finite & (indices >= bin_count)))
    accepted = finite & (indices >= 0) & (indices < bin_count)

    shape = (len(event_weights), bin_count)
    count_matrix = np.zeros(shape, dtype=float)
    response_matrix = np.zeros(shape, dtype=float)
    deficit_matrix = np.zeros(shape, dtype=float)
    pt_matrix = np.zeros(shape, dtype=float)
    accepted_events = event_indices[accepted]
    accepted_bins = indices[accepted]
    np.add.at(count_matrix, (accepted_events, accepted_bins), 1.0)
    np.add.at(
        response_matrix,
        (accepted_events, accepted_bins),
        fractional_shift[accepted],
    )
    np.add.at(
        deficit_matrix,
        (accepted_events, accepted_bins),
        baseline_pt[accepted] * fractional_shift[accepted],
    )
    np.add.at(pt_matrix, (accepted_events, accepted_bins), baseline_pt[accepted])

    weights = np.asarray(event_weights, dtype=float)[:, np.newaxis]
    weighted_counts = count_matrix * weights
    weighted_responses = response_matrix * weights
    weighted_deficits = deficit_matrix * weights
    weighted_pt = pt_matrix * weights
    counts = np.sum(weighted_counts, axis=0)
    responses = np.sum(weighted_responses, axis=0)
    deficits = np.sum(weighted_deficits, axis=0)
    pt_totals = np.sum(weighted_pt, axis=0)
    mean = np.divide(
        responses,
        counts,
        out=np.full(bin_count, np.nan),
        where=counts != 0.0,
    )
    pt_fraction = np.divide(
        deficits,
        pt_totals,
        out=np.full(bin_count, np.nan),
        where=pt_totals != 0.0,
    )

    leave_counts = counts[np.newaxis, :] - weighted_counts
    leave_responses = responses[np.newaxis, :] - weighted_responses
    leave_pt = pt_totals[np.newaxis, :] - weighted_pt
    leave_deficits = deficits[np.newaxis, :] - weighted_deficits
    leave_mean = np.divide(
        leave_responses,
        leave_counts,
        out=np.full_like(leave_responses, np.nan),
        where=leave_counts != 0.0,
    )
    leave_pt_fraction = np.divide(
        leave_deficits,
        leave_pt,
        out=np.full_like(leave_deficits, np.nan),
        where=leave_pt != 0.0,
    )
    return BinnedResponse(
        mean=mean,
        mean_error=np.asarray(jackknife_error(leave_mean)),
        pt_weighted_fraction=pt_fraction,
        pt_weighted_fraction_error=np.asarray(jackknife_error(leave_pt_fraction)),
        raw_entries=np.sum(count_matrix, axis=0).astype(np.int64),
        weighted_entries=counts,
        leave_mean=leave_mean,
        leave_pt_weighted_fraction=leave_pt_fraction,
        underflow_entries=underflow,
        overflow_entries=overflow,
    )


def pt_selection_label(pt_min: float, pt_max: float | None) -> str:
    if pt_max is None:
        return rf"$p_T^{{\rm no-pre}}>{pt_min:g}$ GeV"
    return rf"${pt_min:g}<p_T^{{\rm no-pre}}\leq {pt_max:g}$ GeV"


def format_float(value: float) -> str:
    return "nan" if not math.isfinite(value) else f"{value:.12e}"


def plot_radius(
    out_dir: Path,
    prefix: str,
    radius_digit: int,
    pt_min: float,
    pt_max: float | None,
    abs_eta_max: float,
    results: dict[str, dict[str, BinnedResponse]],
) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    specs = proxy_specs()
    figure, axes = plt.subplots(2, 2, figsize=(13.0, 8.4))
    for axis, spec in zip(axes.flat, specs):
        positions = np.arange(len(spec.bin_labels), dtype=float)
        panel_limits: list[float] = []
        for flavor in FLAVORS:
            result = results[spec.key][flavor]
            finite = np.isfinite(result.mean) & np.isfinite(result.mean_error)
            if np.any(finite):
                panel_limits.extend(
                    (100.0 * (result.mean[finite] - result.mean_error[finite])).tolist()
                )
                panel_limits.extend(
                    (100.0 * (result.mean[finite] + result.mean_error[finite])).tolist()
                )
            axis.errorbar(
                positions[finite],
                100.0 * result.mean[finite],
                yerr=100.0 * result.mean_error[finite],
                color=FLAVOR_COLORS[flavor],
                marker=FLAVOR_MARKERS[flavor],
                markersize=4.0,
                linewidth=1.0,
                capsize=2.0,
                label=FLAVOR_LABELS[flavor],
            )
        axis.axhline(0.0, color="0.45", linewidth=0.9)
        axis.set_xticks(positions)
        axis.set_xticklabels(spec.bin_labels, rotation=30, ha="right")
        axis.set_xlabel(spec.label)
        axis.grid(alpha=0.22)
        if panel_limits:
            low = min(-1.0, min(panel_limits))
            high = max(1.0, max(panel_limits))
            margin = 0.12 * (high - low)
            axis.set_ylim(low - margin, high + margin)
    axes[0, 0].set_ylabel(r"paired early-stage shift $\langle1-p_T^{\rm pre}/p_T^{\rm no}\rangle$ [\%]")
    axes[1, 0].set_ylabel(r"paired early-stage shift $\langle1-p_T^{\rm pre}/p_T^{\rm no}\rangle$ [\%]")
    axes[0, 0].legend(frameon=False, fontsize=9)
    figure.suptitle(
        rf"O16+O16 5.36 TeV, 0--5% diagnostic; anti-$k_T$ $R={radius_digit / 10:.1f}$, "
        + pt_selection_label(pt_min, pt_max)
        + rf", $|\eta^{{\rm no-pre}}|<{abs_eta_max:g}$",
        fontsize=13,
    )
    figure.text(
        0.5,
        0.012,
        "Selection and proxy bins use the no-prehydro jet; matched Plan-B jet supplies the numerator. "
        "Positive values mean additional loss. Event-weighted paired delete-one-run jackknife.",
        ha="center",
        fontsize=8.5,
    )
    figure.subplots_adjust(top=0.92, bottom=0.16, wspace=0.17, hspace=0.39)
    stem = out_dir / f"{prefix}_R0{radius_digit}_charge_response"
    pdf = stem.with_suffix(".pdf")
    png = stem.with_suffix(".png")
    figure.savefig(pdf)
    figure.savefig(png, dpi=180)
    plt.close(figure)
    return pdf, png


def plot_single_many_summary(
    out_dir: Path,
    prefix: str,
    pt_min: float,
    pt_max: float | None,
    abs_eta_max: float,
    category_results: dict[int, dict[str, BinnedResponse]],
) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(1, 2, figsize=(12.2, 4.7))
    positions = np.arange(3, dtype=float)
    for radius_digit in RADIUS_DIGITS:
        result = category_results[radius_digit]["quark"]
        axes[0].errorbar(
            positions,
            100.0 * result.mean,
            yerr=100.0 * result.mean_error,
            color=RADIUS_COLORS[radius_digit],
            marker=RADIUS_MARKERS[radius_digit],
            linewidth=1.2,
            capsize=2.0,
            label=f"R=0.{radius_digit}",
        )
    axes[0].axhline(0.0, color="0.45", linewidth=0.9)
    axes[0].set_xticks(positions)
    axes[0].set_xticklabels(CATEGORY_LABELS)
    axes[0].set_ylabel(r"quark-tagged $\langle1-p_T^{\rm pre}/p_T^{\rm no}\rangle$ [\%]")
    axes[0].set_title("Wake-excluded effective-multiplicity classes")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.22)

    radius_positions = np.arange(len(RADIUS_DIGITS), dtype=float)
    for flavor in FLAVORS:
        differences = []
        errors = []
        for radius_digit in RADIUS_DIGITS:
            result = category_results[radius_digit][flavor]
            differences.append(result.mean[0] - result.mean[2])
            errors.append(
                jackknife_error(result.leave_mean[:, 0] - result.leave_mean[:, 2])
            )
        axes[1].errorbar(
            radius_positions,
            100.0 * np.asarray(differences),
            yerr=100.0 * np.asarray(errors),
            color=FLAVOR_COLORS[flavor],
            marker=FLAVOR_MARKERS[flavor],
            linewidth=1.2,
            capsize=2.0,
            label=FLAVOR_LABELS[flavor],
        )
    axes[1].axhline(0.0, color="0.45", linewidth=0.9)
    axes[1].set_xticks(radius_positions)
    axes[1].set_xticklabels([f"R=0.{value}" for value in RADIUS_DIGITS])
    axes[1].set_ylabel("single-like minus many-like shift [%]")
    axes[1].set_title("Positive values support the proposed ordering")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.22)

    figure.suptitle(
        "O16+O16 5.36 TeV, 0--5% diagnostic; "
        + pt_selection_label(pt_min, pt_max)
        + rf", $|\eta^{{\rm no-pre}}|<{abs_eta_max:g}$",
        fontsize=13,
    )
    figure.text(
        0.5,
        0.012,
        "No-prehydro categories; matched Plan-B response. Quark/gluon subsets require the same hard-parton marker in both variants.",
        ha="center",
        fontsize=8.5,
    )
    figure.subplots_adjust(top=0.84, bottom=0.19, wspace=0.27)
    stem = out_dir / f"{prefix}_quark_single_many_summary"
    pdf = stem.with_suffix(".pdf")
    png = stem.with_suffix(".png")
    figure.savefig(pdf)
    figure.savefig(png, dpi=180)
    plt.close(figure)
    return pdf, png


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--pt-min", type=float, default=30.0)
    parser.add_argument("--pt-max", type=float)
    parser.add_argument("--abs-eta-max", type=float, default=2.0)
    parser.add_argument(
        "--max-pair-match-dr-fraction",
        type=float,
        default=0.5,
        help="require paired jet axes to satisfy DeltaR < fraction*R",
    )
    parser.add_argument("--prefix", default="oo5360_jet_charge_response_pt30")
    return parser


def run(args: argparse.Namespace) -> dict[str, object]:
    if args.pt_min <= 0.0:
        raise ValueError("--pt-min must be positive")
    if args.pt_max is not None and args.pt_max <= args.pt_min:
        raise ValueError("--pt-max must be greater than --pt-min")
    if args.abs_eta_max <= 0.0:
        raise ValueError("--abs-eta-max must be positive")
    if not 0.0 < args.max_pair_match_dr_fraction <= 0.5:
        raise ValueError("--max-pair-match-dr-fraction must be in (0, 0.5]")
    input_root = args.input_root.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    with uproot.open(input_root) as root_file:
        pairs = root_file["Pairs"].arrays(
            ["pairId", "eventWeight", "seed"], library="np"
        )
        pair_ids = pairs["pairId"]
        event_weights = pairs["eventWeight"].astype(float)
        if len(pair_ids) == 0 or len(np.unique(pair_ids)) != len(pair_ids):
            raise ValueError("Pairs tree has empty or duplicate pairId values")
        if len(np.unique(pairs["seed"])) != len(pair_ids):
            raise ValueError("Pairs tree contains duplicate hard-event seeds")
        if not np.all(np.isfinite(event_weights)) or not np.all(event_weights > 0.0):
            raise ValueError("event weights must be finite and positive")

        branches = ["pairId", "eventWeight"]
        for radius_digit in RADIUS_DIGITS:
            branches.extend(
                f"jet{radius_digit}{suffix}"
                for suffix in (
                    "Pt",
                    "Eta",
                    "NormalEffectiveMultiplicity",
                    "LeadingNormalFraction",
                    "NSD",
                    "MaxKt",
                    "HardPartonId",
                    "PairMatchIndex",
                    "PairMatchDR",
                    "PairMatchOtherPt",
                    "PairMatchOtherHardPartonId",
                )
            )
        arrays = root_file["noPrehydro/Jets"].arrays(branches, library="ak")
    if not np.array_equal(ak.to_numpy(arrays.pairId), pair_ids):
        raise ValueError("noPrehydro/Jets pairId ordering differs from Pairs")
    if not np.array_equal(ak.to_numpy(arrays.eventWeight), event_weights):
        raise ValueError("noPrehydro/Jets eventWeight differs from Pairs")

    histogram_rows: list[list[object]] = []
    summary_rows: list[list[object]] = []
    metadata: dict[str, object] = {
        "inputRoot": str(input_root),
        "inputBytes": input_root.stat().st_size,
        "pairCount": len(pair_ids),
        "ptMinGeV": args.pt_min,
        "ptMaxGeV": args.pt_max,
        "absEtaMax": args.abs_eta_max,
        "maxPairMatchDRFraction": args.max_pair_match_dr_fraction,
        "selection": "no-prehydro corrected jet pT and eta; paired-axis match required",
        "fractionalShift": "1 - matched Plan-B corrected pT / no-prehydro corrected pT",
        "positiveShiftMeaning": "additional loss when Plan-B prehydro is enabled",
        "normalEffectiveMultiplicity": "(sum normal-hadron pT)^2 / sum normal-hadron pT^2; wake excluded",
        "hardPartonTag": "nearest one-to-one raw-label -2 outgoing marker within DeltaR<R; flavor subsets require the same marker PDG ID in both variants",
        "normalization": "conditional means use the biased-PYTHIA event weight once; common merged cross-section factor cancels",
        "uncertainty": "paired delete-one-run jackknife",
        "truthCaveat": "final-state proxies do not measure the number of active shower charges at hydro start",
        "radii": {},
    }
    category_results: dict[int, dict[str, BinnedResponse]] = {}

    for radius_digit in RADIUS_DIGITS:
        prefix = f"jet{radius_digit}"
        pt = arrays[f"{prefix}Pt"]
        eta = arrays[f"{prefix}Eta"]
        match_index = arrays[f"{prefix}PairMatchIndex"]
        match_dr = arrays[f"{prefix}PairMatchDR"]
        other_pt = arrays[f"{prefix}PairMatchOtherPt"]
        selected = (
            np.isfinite(pt)
            & (pt > args.pt_min)
            & np.isfinite(eta)
            & (abs(eta) < args.abs_eta_max)
            & (match_index >= 0)
            & np.isfinite(match_dr)
            & (match_dr < args.max_pair_match_dr_fraction * radius_digit / 10.0)
            & np.isfinite(other_pt)
            & (other_pt > 0.0)
        )
        if args.pt_max is not None:
            selected = selected & (pt <= args.pt_max)
        event_lengths = ak.to_numpy(ak.num(pt[selected], axis=1))
        event_indices = np.repeat(np.arange(len(event_weights), dtype=np.int64), event_lengths)
        baseline_pt = ak.to_numpy(ak.flatten(pt[selected])).astype(float)
        matched_pt = ak.to_numpy(ak.flatten(other_pt[selected])).astype(float)
        fractional_shift = 1.0 - matched_pt / baseline_pt
        hard_ids = ak.to_numpy(ak.flatten(arrays[f"{prefix}HardPartonId"][selected])).astype(int)
        other_hard_ids = ak.to_numpy(
            ak.flatten(arrays[f"{prefix}PairMatchOtherHardPartonId"][selected])
        ).astype(int)
        proxy_values = {
            "normal_neff": ak.to_numpy(
                ak.flatten(arrays[f"{prefix}NormalEffectiveMultiplicity"][selected])
            ).astype(float),
            "leading_normal_fraction": ak.to_numpy(
                ak.flatten(arrays[f"{prefix}LeadingNormalFraction"][selected])
            ).astype(float),
            "resolved_primary_prongs": 1.0
            + ak.to_numpy(ak.flatten(arrays[f"{prefix}NSD"][selected])).astype(float),
            "max_kt": np.nan_to_num(
                ak.to_numpy(ak.flatten(arrays[f"{prefix}MaxKt"][selected])).astype(float),
                nan=0.0,
            ),
        }

        tagged = hard_ids != 0
        shared_tag = tagged & (other_hard_ids == hard_ids)
        radius_metadata: dict[str, object] = {
            "matchedSelectedJets": len(baseline_pt),
            "hardPartonTaggedJets": int(np.count_nonzero(tagged)),
            "sharedHardPartonTaggedJets": int(np.count_nonzero(shared_tag)),
            "hardPartonTagFraction": (
                float(np.count_nonzero(tagged) / len(tagged)) if len(tagged) else math.nan
            ),
            "quarkTaggedJets": int(
                np.count_nonzero(flavor_mask(hard_ids, "quark", other_hard_ids))
            ),
            "gluonTaggedJets": int(
                np.count_nonzero(flavor_mask(hard_ids, "gluon", other_hard_ids))
            ),
            "observables": {},
        }
        plot_results: dict[str, dict[str, BinnedResponse]] = {}
        for spec in proxy_specs():
            plot_results[spec.key] = {}
            observable_metadata: dict[str, object] = {}
            for flavor in FLAVORS:
                keep = flavor_mask(hard_ids, flavor, other_hard_ids)
                result = binned_response(
                    proxy_values[spec.key][keep],
                    fractional_shift[keep],
                    baseline_pt[keep],
                    event_indices[keep],
                    event_weights,
                    spec.edges,
                )
                plot_results[spec.key][flavor] = result
                observable_metadata[flavor] = {
                    "rawEntriesInRange": int(np.sum(result.raw_entries)),
                    "underflowEntries": result.underflow_entries,
                    "overflowEntries": result.overflow_entries,
                }
                for bin_index in range(len(spec.edges) - 1):
                    histogram_rows.append(
                        [
                            f"0.{radius_digit}",
                            flavor,
                            spec.key,
                            bin_index,
                            format_float(float(spec.edges[bin_index])),
                            format_float(float(spec.edges[bin_index + 1])),
                            spec.bin_labels[bin_index],
                            int(result.raw_entries[bin_index]),
                            format_float(float(result.weighted_entries[bin_index])),
                            format_float(float(result.mean[bin_index])),
                            format_float(float(result.mean_error[bin_index])),
                            format_float(float(result.pt_weighted_fraction[bin_index])),
                            format_float(float(result.pt_weighted_fraction_error[bin_index])),
                        ]
                    )
            radius_metadata["observables"][spec.key] = observable_metadata

        category_edges = np.array([0.5, 3.0, 8.0, np.inf])
        category_labels = ("single_like_neff_lt3", "intermediate_neff_3to8", "many_like_neff_ge8")
        category_metadata: dict[str, object] = {}
        category_results[radius_digit] = {}
        for flavor in FLAVORS:
            keep = flavor_mask(hard_ids, flavor, other_hard_ids)
            category = binned_response(
                proxy_values["normal_neff"][keep],
                fractional_shift[keep],
                baseline_pt[keep],
                event_indices[keep],
                event_weights,
                category_edges,
            )
            category_results[radius_digit][flavor] = category
            difference = float(category.mean[0] - category.mean[2])
            difference_error = float(
                jackknife_error(category.leave_mean[:, 0] - category.leave_mean[:, 2])
            )
            category_metadata[flavor] = {
                "singleMinusManyMeanFractionalShift": difference,
                "singleMinusManyStatError": difference_error,
                "hypothesisExpectsPositive": True,
            }
            for index, label in enumerate(category_labels):
                summary_rows.append(
                    [
                        f"0.{radius_digit}",
                        flavor,
                        label,
                        int(category.raw_entries[index]),
                        format_float(float(category.weighted_entries[index])),
                        format_float(float(category.mean[index])),
                        format_float(float(category.mean_error[index])),
                        format_float(float(category.pt_weighted_fraction[index])),
                        format_float(float(category.pt_weighted_fraction_error[index])),
                        format_float(difference),
                        format_float(difference_error),
                    ]
                )
        radius_metadata["singleManyComparison"] = category_metadata
        pdf, png = plot_radius(
            out_dir,
            args.prefix,
            radius_digit,
            args.pt_min,
            args.pt_max,
            args.abs_eta_max,
            plot_results,
        )
        radius_metadata["plotPdf"] = str(pdf)
        radius_metadata["plotPng"] = str(png)
        metadata["radii"][f"0.{radius_digit}"] = radius_metadata

    summary_pdf, summary_png = plot_single_many_summary(
        out_dir,
        args.prefix,
        args.pt_min,
        args.pt_max,
        args.abs_eta_max,
        category_results,
    )
    metadata["singleManySummaryPlotPdf"] = str(summary_pdf)
    metadata["singleManySummaryPlotPng"] = str(summary_png)

    histogram_path = out_dir / f"{args.prefix}_charge_response.tsv"
    with histogram_path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "radius",
                "flavor",
                "observable",
                "bin_index",
                "bin_low",
                "bin_high",
                "bin_label",
                "raw_jets",
                "weighted_jets",
                "mean_fractional_early_shift",
                "mean_stat_error",
                "pt_weighted_fractional_early_shift",
                "pt_weighted_stat_error",
            ]
        )
        writer.writerows(histogram_rows)

    summary_path = out_dir / f"{args.prefix}_single_many_summary.tsv"
    with summary_path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "radius",
                "flavor",
                "category",
                "raw_jets",
                "weighted_jets",
                "mean_fractional_early_shift",
                "mean_stat_error",
                "pt_weighted_fractional_early_shift",
                "pt_weighted_stat_error",
                "single_minus_many_mean_shift",
                "single_minus_many_stat_error",
            ]
        )
        writer.writerows(summary_rows)

    metadata["histogramTsv"] = str(histogram_path)
    metadata["singleManySummaryTsv"] = str(summary_path)
    metadata_path = out_dir / f"{args.prefix}_charge_response_metadata.json"
    metadata["metadataJson"] = str(metadata_path)
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return metadata


def main() -> int:
    args = make_parser().parse_args()
    metadata = run(args)
    print(json.dumps(metadata, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
