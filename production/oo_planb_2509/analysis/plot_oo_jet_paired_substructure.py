#!/usr/bin/env python3
"""Measure migration-safe paired OO jet substructure changes."""

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


RADIUS_DIGITS = (1, 2, 4, 8)
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
TRANSITIONS = (
    ("fail_to_fail", False, False),
    ("fail_to_pass", False, True),
    ("pass_to_fail", True, False),
    ("pass_to_pass", True, True),
)
CATEGORY_LABELS = {
    "all": "All",
    "no_sd_fail": "SD fail",
    "no_sd_pass": "SD pass",
    "single_like": "Single-like",
    "intermediate": "Intermediate",
    "many_like": "Many-like",
}
PLOT_CATEGORIES = tuple(CATEGORY_LABELS)
OBSERVABLES = (
    "epsilon_pre",
    "pt_weighted_epsilon_pre",
    "delta_ptd",
    "relative_ptd",
    "relative_mean_ptd_shift",
    "delta_normal_ptd",
    "relative_normal_ptd",
    "relative_mean_normal_ptd_shift",
    "delta_mult",
    "relative_mean_mult_shift",
    "delta_total_mult",
    "relative_mean_total_mult_shift",
    "delta_zg",
    "relative_mean_zg_shift",
    "delta_rg",
    "relative_mean_rg_shift",
    "delta_girth",
    "relative_mean_girth_shift",
    "delta_maxkt",
    "relative_mean_maxkt_shift",
    "delta_mean_tauf",
    "relative_mean_tauf_shift",
    "delta_mean_log10_tauf",
    "delta_hardest_log10_tauf",
)


@dataclasses.dataclass
class RatioEstimate:
    value: float
    event_error: float
    hydro_error: float
    raw_entries: int
    raw_events: int
    weighted_entries: float
    normalization_denominator: float
    effective_event_contributors: float
    max_event_fraction: float
    leave_event: np.ndarray
    leave_hydro: np.ndarray


@dataclasses.dataclass
class DifferenceEstimate:
    value: float
    event_error: float
    hydro_error: float


def format_float(value: float) -> str:
    return "nan" if not math.isfinite(value) else f"{value:.12e}"


def jackknife_error(leave_one_out: np.ndarray) -> float:
    values = np.asarray(leave_one_out, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.all(np.isfinite(values)):
        return math.nan
    mean = float(np.mean(values))
    return math.sqrt(
        (len(values) - 1.0) / len(values)
        * float(np.sum((values - mean) ** 2))
    )


def ratio_estimate(
    numerator: np.ndarray,
    denominator: np.ndarray,
    entry_mask: np.ndarray,
    event_indices: np.ndarray,
    event_weights: np.ndarray,
    event_hydro_indices: np.ndarray,
) -> RatioEstimate:
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    entry_mask = np.asarray(entry_mask, dtype=bool)
    event_indices = np.asarray(event_indices, dtype=np.int64)
    if not (
        len(numerator)
        == len(denominator)
        == len(entry_mask)
        == len(event_indices)
    ):
        raise ValueError("jet contribution arrays have inconsistent lengths")
    if np.any(event_indices < 0) or np.any(event_indices >= len(event_weights)):
        raise ValueError("jet event index is outside the Pairs tree")

    finite = np.isfinite(numerator) & np.isfinite(denominator) & entry_mask
    numerator = np.where(finite, numerator, 0.0)
    denominator = np.where(finite, denominator, 0.0)
    entries = finite.astype(float)

    event_numerator = np.zeros(len(event_weights), dtype=float)
    event_denominator = np.zeros(len(event_weights), dtype=float)
    event_entries = np.zeros(len(event_weights), dtype=float)
    np.add.at(event_numerator, event_indices, numerator)
    np.add.at(event_denominator, event_indices, denominator)
    np.add.at(event_entries, event_indices, entries)

    weighted_numerator = event_weights * event_numerator
    weighted_denominator = event_weights * event_denominator
    weighted_entries_by_event = event_weights * event_entries
    numerator_total = float(np.sum(weighted_numerator))
    denominator_total = float(np.sum(weighted_denominator))
    value = (
        numerator_total / denominator_total
        if denominator_total != 0.0
        else math.nan
    )

    leave_denominator = denominator_total - weighted_denominator
    leave_event = np.divide(
        numerator_total - weighted_numerator,
        leave_denominator,
        out=np.full(len(event_weights), np.nan),
        where=leave_denominator != 0.0,
    )

    hydro_values, event_to_hydro = np.unique(
        event_hydro_indices, return_inverse=True
    )
    hydro_numerator = np.zeros(len(hydro_values), dtype=float)
    hydro_denominator = np.zeros(len(hydro_values), dtype=float)
    np.add.at(hydro_numerator, event_to_hydro, weighted_numerator)
    np.add.at(hydro_denominator, event_to_hydro, weighted_denominator)
    leave_hydro_denominator = denominator_total - hydro_denominator
    leave_hydro = np.divide(
        numerator_total - hydro_numerator,
        leave_hydro_denominator,
        out=np.full(len(hydro_values), np.nan),
        where=leave_hydro_denominator != 0.0,
    )

    weighted_entries = float(np.sum(weighted_entries_by_event))
    square_sum = float(
        np.sum(weighted_entries_by_event * weighted_entries_by_event)
    )
    effective = (
        weighted_entries * weighted_entries / square_sum
        if square_sum > 0.0
        else 0.0
    )
    max_fraction = (
        float(np.max(weighted_entries_by_event) / weighted_entries)
        if weighted_entries > 0.0
        else math.nan
    )
    return RatioEstimate(
        value=value,
        event_error=jackknife_error(leave_event),
        hydro_error=jackknife_error(leave_hydro),
        raw_entries=int(np.sum(entries)),
        raw_events=int(np.count_nonzero(event_entries)),
        weighted_entries=weighted_entries,
        normalization_denominator=denominator_total,
        effective_event_contributors=effective,
        max_event_fraction=max_fraction,
        leave_event=leave_event,
        leave_hydro=leave_hydro,
    )


def fraction_estimate(
    numerator_mask: np.ndarray,
    denominator_mask: np.ndarray,
    event_indices: np.ndarray,
    event_weights: np.ndarray,
    event_hydro_indices: np.ndarray,
) -> RatioEstimate:
    return ratio_estimate(
        numerator=np.asarray(numerator_mask, dtype=float),
        denominator=np.asarray(denominator_mask, dtype=float),
        entry_mask=np.asarray(denominator_mask, dtype=bool),
        event_indices=event_indices,
        event_weights=event_weights,
        event_hydro_indices=event_hydro_indices,
    )


def mean_estimate(
    values: np.ndarray,
    mask: np.ndarray,
    event_indices: np.ndarray,
    event_weights: np.ndarray,
    event_hydro_indices: np.ndarray,
) -> RatioEstimate:
    values = np.asarray(values, dtype=float)
    mask = np.asarray(mask, dtype=bool) & np.isfinite(values)
    return ratio_estimate(
        numerator=values,
        denominator=np.ones(len(values), dtype=float),
        entry_mask=mask,
        event_indices=event_indices,
        event_weights=event_weights,
        event_hydro_indices=event_hydro_indices,
    )


def difference_estimate(
    minuend: RatioEstimate, subtrahend: RatioEstimate
) -> DifferenceEstimate:
    if len(minuend.leave_event) != len(subtrahend.leave_event):
        raise ValueError("event jackknife arrays have inconsistent lengths")
    if len(minuend.leave_hydro) != len(subtrahend.leave_hydro):
        raise ValueError("hydro jackknife arrays have inconsistent lengths")
    return DifferenceEstimate(
        value=minuend.value - subtrahend.value,
        event_error=jackknife_error(
            minuend.leave_event - subtrahend.leave_event
        ),
        hydro_error=jackknife_error(
            minuend.leave_hydro - subtrahend.leave_hydro
        ),
    )


def flavor_mask(
    hard_parton_id: np.ndarray,
    matched_hard_parton_id: np.ndarray,
    flavor: str,
) -> np.ndarray:
    if flavor == "all":
        return np.ones(len(hard_parton_id), dtype=bool)
    shared = (hard_parton_id != 0) & (
        hard_parton_id == matched_hard_parton_id
    )
    if flavor == "quark":
        absolute = np.abs(hard_parton_id)
        return shared & (absolute >= 1) & (absolute <= 6)
    if flavor == "gluon":
        return shared & (hard_parton_id == 21)
    raise ValueError(f"unsupported flavor {flavor}")


def category_masks(
    no_soft_drop_valid: np.ndarray,
    normal_effective_multiplicity: np.ndarray,
) -> dict[str, np.ndarray]:
    valid_neff = np.isfinite(normal_effective_multiplicity)
    return {
        "all": np.ones(len(no_soft_drop_valid), dtype=bool),
        "no_sd_fail": ~no_soft_drop_valid,
        "no_sd_pass": no_soft_drop_valid,
        "single_like": valid_neff & (normal_effective_multiplicity < 3.0),
        "intermediate": (
            valid_neff
            & (normal_effective_multiplicity >= 3.0)
            & (normal_effective_multiplicity < 8.0)
        ),
        "many_like": valid_neff & (normal_effective_multiplicity >= 8.0),
    }


def selected_flat(jagged_values: ak.Array, selection: ak.Array) -> np.ndarray:
    return ak.to_numpy(ak.flatten(jagged_values[selection]))


def matched_flat(
    jagged_values: ak.Array,
    safe_match_index: ak.Array,
    selection: ak.Array,
) -> np.ndarray:
    return ak.to_numpy(
        ak.flatten(jagged_values[safe_match_index[selection]])
    )


def formation_jet_summaries(
    flat_tau_f: ak.Array, offsets: ak.Array
) -> tuple[ak.Array, ak.Array]:
    """Return per-jet arithmetic mean tau_f and mean log10(tau_f)."""
    mean_tau_events: list[list[float]] = []
    mean_log_events: list[list[float]] = []
    for event_values, event_offsets in zip(
        ak.to_list(flat_tau_f), ak.to_list(offsets), strict=True
    ):
        if not event_offsets or event_offsets[0] != 0:
            raise ValueError("formation offsets must start at zero")
        if event_offsets[-1] != len(event_values):
            raise ValueError("formation offsets do not close on the flat split vector")
        if any(second < first for first, second in zip(event_offsets, event_offsets[1:])):
            raise ValueError("formation offsets must be monotonic")
        event_means: list[float] = []
        event_log_means: list[float] = []
        for first, second in zip(event_offsets, event_offsets[1:]):
            values = np.asarray(event_values[first:second], dtype=float)
            values = values[np.isfinite(values) & (values > 0.0)]
            if values.size:
                event_means.append(float(np.mean(values)))
                event_log_means.append(float(np.mean(np.log10(values))))
            else:
                event_means.append(math.nan)
                event_log_means.append(math.nan)
        mean_tau_events.append(event_means)
        mean_log_events.append(event_log_means)
    return ak.Array(mean_tau_events), ak.Array(mean_log_events)


def estimate_to_dict(estimate: RatioEstimate) -> dict[str, object]:
    return {
        "value": estimate.value,
        "eventJackknifeError": estimate.event_error,
        "hydroBlockJackknifeError": estimate.hydro_error,
        "rawEntries": estimate.raw_entries,
        "rawEvents": estimate.raw_events,
        "weightedEntries": estimate.weighted_entries,
        "normalizationDenominator": estimate.normalization_denominator,
        "effectiveEventContributors": estimate.effective_event_contributors,
        "maxEventContributionFraction": estimate.max_event_fraction,
    }


def difference_to_dict(estimate: DifferenceEstimate) -> dict[str, float]:
    return {
        "value": estimate.value,
        "eventJackknifeError": estimate.event_error,
        "hydroBlockJackknifeError": estimate.hydro_error,
        "eventSignificance": (
            estimate.value / estimate.event_error
            if estimate.event_error > 0.0
            else math.nan
        ),
        "hydroSignificance": (
            estimate.value / estimate.hydro_error
            if estimate.hydro_error > 0.0
            else math.nan
        ),
    }


def finite_errorbar(axis, x, values, errors, **kwargs) -> None:
    x = np.asarray(x, dtype=float)
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    finite = np.isfinite(values) & np.isfinite(errors)
    axis.errorbar(x[finite], values[finite], yerr=errors[finite], **kwargs)


def draw_transition_matrix(axis, estimates, title: str) -> None:
    values = np.array(
        [
            [estimates["fail_to_fail"].value, estimates["fail_to_pass"].value],
            [estimates["pass_to_fail"].value, estimates["pass_to_pass"].value],
        ]
    )
    errors = np.array(
        [
            [
                estimates["fail_to_fail"].event_error,
                estimates["fail_to_pass"].event_error,
            ],
            [
                estimates["pass_to_fail"].event_error,
                estimates["pass_to_pass"].event_error,
            ],
        ]
    )
    axis.imshow(
        np.log10(np.clip(values, 1.0e-6, None)),
        cmap="Blues",
        vmin=-4.0,
        vmax=0.0,
    )
    for row in range(2):
        for column in range(2):
            value = values[row, column]
            error = errors[row, column]
            color = "white" if value > 0.2 else "#222222"
            axis.text(
                column,
                row,
                f"{100.0 * value:.3f}%\n"
                f"+/- {100.0 * error:.3f}%",
                ha="center",
                va="center",
                fontsize=8,
                color=color,
            )
    axis.set_xticks((0, 1), ("Pre-hydro fail", "Pre-hydro pass"))
    axis.set_yticks((0, 1), ("No-pre fail", "No-pre pass"))
    axis.set_title(title, fontsize=10)
    axis.tick_params(labelsize=8)


def plot_radius(
    out_dir: Path,
    prefix: str,
    radius_digit: int,
    pt_min: float,
    pt_max: float | None,
    transitions: dict[str, dict[str, RatioEstimate]],
    fail_summary: dict[str, dict[str, RatioEstimate]],
    observables: dict[str, dict[str, dict[str, RatioEstimate]]],
) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 3, figsize=(15.4, 8.3))
    draw_transition_matrix(
        axes[0, 0], transitions["all"], "Soft Drop transitions: all jets"
    )
    draw_transition_matrix(
        axes[0, 1], transitions["quark"], "Soft Drop transitions: quark tagged"
    )

    x_flavor = np.arange(len(FLAVORS), dtype=float)
    for offset, key, label, marker in (
        (-0.08, "no_fail_fraction", "No-pre fail", "o"),
        (0.08, "pre_fail_fraction", "Pre-hydro fail", "s"),
    ):
        finite_errorbar(
            axes[0, 2],
            x_flavor + offset,
            [100.0 * fail_summary[item][key].value for item in FLAVORS],
            [
                100.0 * fail_summary[item][key].event_error
                for item in FLAVORS
            ],
            marker=marker,
            linestyle="none",
            capsize=2.5,
            label=label,
        )
    axes[0, 2].set_xticks(
        x_flavor, [FLAVOR_LABELS[item] for item in FLAVORS]
    )
    axes[0, 2].set_ylabel("failed-Soft-Drop fraction [%]")
    axes[0, 2].legend(frameon=False, fontsize=8)
    axes[0, 2].grid(alpha=0.22)

    finite_errorbar(
        axes[1, 0],
        x_flavor,
        [
            100.0 * fail_summary[item]["delta_fail_fraction"].value
            for item in FLAVORS
        ],
        [
            100.0 * fail_summary[item]["delta_fail_fraction"].event_error
            for item in FLAVORS
        ],
        marker="o",
        linestyle="none",
        color="#333333",
        capsize=2.5,
    )
    axes[1, 0].axhline(0.0, color="0.5", linewidth=0.9)
    axes[1, 0].set_xticks(
        x_flavor, [FLAVOR_LABELS[item] for item in FLAVORS]
    )
    axes[1, 0].set_ylabel("Pre-hydro minus no-pre fail fraction [points]")
    axes[1, 0].grid(alpha=0.22)

    x_category = np.arange(len(PLOT_CATEGORIES), dtype=float)
    quark = observables["quark"]
    finite_errorbar(
        axes[1, 1],
        x_category,
        [
            100.0 * quark[category]["epsilon_pre"].value
            for category in PLOT_CATEGORIES
        ],
        [
            100.0 * quark[category]["epsilon_pre"].event_error
            for category in PLOT_CATEGORIES
        ],
        marker="o",
        linestyle="none",
        color=FLAVOR_COLORS["quark"],
        capsize=2.5,
    )
    axes[1, 1].axhline(0.0, color="0.5", linewidth=0.9)
    axes[1, 1].set_xticks(
        x_category,
        [CATEGORY_LABELS[item] for item in PLOT_CATEGORIES],
        rotation=25,
        ha="right",
    )
    axes[1, 1].set_ylabel(r"quark-tagged $\epsilon_{\rm pre}$ [%]")
    axes[1, 1].grid(alpha=0.22)

    for offset, key, label, marker, color in (
        (-0.08, "delta_ptd", r"$\Delta p_{TD}$", "o", "#0072B2"),
        (
            0.08,
            "delta_normal_ptd",
            r"$\Delta p_{TD}^{\rm normal}$",
            "s",
            "#D55E00",
        ),
    ):
        finite_errorbar(
            axes[1, 2],
            x_category + offset,
            [quark[category][key].value for category in PLOT_CATEGORIES],
            [quark[category][key].event_error for category in PLOT_CATEGORIES],
            marker=marker,
            linestyle="none",
            color=color,
            capsize=2.5,
            label=label,
        )
    axes[1, 2].axhline(0.0, color="0.5", linewidth=0.9)
    axes[1, 2].set_xticks(
        x_category,
        [CATEGORY_LABELS[item] for item in PLOT_CATEGORIES],
        rotation=25,
        ha="right",
    )
    axes[1, 2].set_ylabel("paired momentum-dispersion shift")
    axes[1, 2].legend(frameon=False, fontsize=8)
    axes[1, 2].grid(alpha=0.22)

    radius = radius_digit / 10.0
    if pt_max is None:
        pt_label = rf"$p_T^{{\rm no}}>{pt_min:g}$ GeV"
    else:
        pt_label = (
            f"{pt_min:g}<"
            + r"$p_T^{\rm no}\leq$"
            + f"{pt_max:g} GeV"
        )
    figure.suptitle(
        rf"O16+O16 5.36 TeV, anti-$k_T$ R={radius:.1f}, {pt_label}; "
        r"select no-pre jet, follow matched pre-hydro jet",
        fontsize=13,
        y=0.99,
    )
    figure.text(
        0.5,
        0.008,
        "Biased-PYTHIA event weight applied once; paired delete-one-run "
        "jackknife shown. Hydro-block errors are stored in TSV/JSON.",
        ha="center",
        fontsize=8.5,
    )
    figure.subplots_adjust(
        left=0.075,
        right=0.98,
        bottom=0.14,
        top=0.91,
        wspace=0.31,
        hspace=0.34,
    )
    stem = out_dir / f"{prefix}_R0{radius_digit}_paired_substructure"
    pdf = stem.with_suffix(".pdf")
    png = stem.with_suffix(".png")
    figure.savefig(pdf)
    figure.savefig(png, dpi=180)
    plt.close(figure)
    return pdf, png


def plot_radius_summary(
    out_dir: Path,
    prefix: str,
    fail_by_radius,
    observable_by_radius,
) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.arange(len(RADIUS_DIGITS), dtype=float)
    labels = [f"R=0.{item}" for item in RADIUS_DIGITS]
    figure, axes = plt.subplots(1, 3, figsize=(14.8, 4.6))

    for offset, flavor in zip((-0.08, 0.08), ("all", "quark")):
        finite_errorbar(
            axes[0],
            x + offset,
            [
                100.0
                * fail_by_radius[radius][flavor][
                    "delta_fail_fraction"
                ].value
                for radius in RADIUS_DIGITS
            ],
            [
                100.0
                * fail_by_radius[radius][flavor][
                    "delta_fail_fraction"
                ].event_error
                for radius in RADIUS_DIGITS
            ],
            marker="o" if flavor == "all" else "s",
            linestyle="none",
            color=FLAVOR_COLORS[flavor],
            capsize=2.5,
            label=FLAVOR_LABELS[flavor],
        )
    axes[0].axhline(0.0, color="0.5", linewidth=0.9)
    axes[0].set_xticks(x, labels)
    axes[0].set_ylabel("Pre-hydro minus no-pre SD-fail fraction [points]")
    axes[0].legend(frameon=False, fontsize=8)
    axes[0].grid(alpha=0.22)

    quark_all = [
        observable_by_radius[radius]["quark"]["all"]
        for radius in RADIUS_DIGITS
    ]
    for offset, key, label, marker, color in (
        (-0.08, "delta_ptd", r"$\Delta p_{TD}$", "o", "#0072B2"),
        (
            0.08,
            "delta_normal_ptd",
            r"$\Delta p_{TD}^{\rm normal}$",
            "s",
            "#D55E00",
        ),
    ):
        finite_errorbar(
            axes[1],
            x + offset,
            [item[key].value for item in quark_all],
            [item[key].event_error for item in quark_all],
            marker=marker,
            linestyle="none",
            color=color,
            capsize=2.5,
            label=label,
        )
    axes[1].axhline(0.0, color="0.5", linewidth=0.9)
    axes[1].set_xticks(x, labels)
    axes[1].set_ylabel("quark-tagged paired dispersion shift")
    axes[1].legend(frameon=False, fontsize=8)
    axes[1].grid(alpha=0.22)

    for offset, category, label, marker, color in (
        (-0.08, "no_sd_fail", "No-pre SD fail", "o", "#CC79A7"),
        (0.08, "no_sd_pass", "No-pre SD pass", "s", "#009E73"),
    ):
        finite_errorbar(
            axes[2],
            x + offset,
            [
                100.0
                * observable_by_radius[radius]["quark"][category][
                    "epsilon_pre"
                ].value
                for radius in RADIUS_DIGITS
            ],
            [
                100.0
                * observable_by_radius[radius]["quark"][category][
                    "epsilon_pre"
                ].event_error
                for radius in RADIUS_DIGITS
            ],
            marker=marker,
            linestyle="none",
            color=color,
            capsize=2.5,
            label=label,
        )
    axes[2].axhline(0.0, color="0.5", linewidth=0.9)
    axes[2].set_xticks(x, labels)
    axes[2].set_ylabel(r"quark-tagged $\epsilon_{\rm pre}$ [%]")
    axes[2].legend(frameon=False, fontsize=8)
    axes[2].grid(alpha=0.22)

    figure.suptitle(
        "Migration-safe paired Soft Drop and momentum-dispersion summary",
        fontsize=13,
    )
    figure.text(
        0.5,
        0.01,
        "Selection and category use the no-prehydro jet; event-level paired "
        "jackknife errors.",
        ha="center",
        fontsize=8.5,
    )
    figure.subplots_adjust(
        left=0.07, right=0.985, bottom=0.16, top=0.88, wspace=0.32
    )
    stem = out_dir / f"{prefix}_paired_substructure_summary"
    pdf = stem.with_suffix(".pdf")
    png = stem.with_suffix(".png")
    figure.savefig(pdf)
    figure.savefig(png, dpi=180)
    plt.close(figure)
    return pdf, png


def plot_matched_observable_summary(
    out_dir: Path,
    prefix: str,
    observable_by_radius,
) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.arange(len(RADIUS_DIGITS), dtype=float)
    labels = [f"R=0.{item}" for item in RADIUS_DIGITS]
    all_jets = {
        radius: observable_by_radius[radius]["all"]["all"]
        for radius in RADIUS_DIGITS
    }
    figure, axes = plt.subplots(2, 3, figsize=(15.4, 8.2))

    for offset, key, label, marker, color in (
        (-0.08, "delta_mult", r"$\Delta N_{+}$", "o", "#0072B2"),
        (0.08, "delta_total_mult", r"$\Delta N_{\rm signed}$", "s", "#D55E00"),
    ):
        finite_errorbar(
            axes[0, 0],
            x + offset,
            [all_jets[radius][key].value for radius in RADIUS_DIGITS],
            [all_jets[radius][key].event_error for radius in RADIUS_DIGITS],
            marker=marker,
            linestyle="none",
            color=color,
            capsize=2.5,
            label=label,
        )
    axes[0, 0].set_ylabel("Pre-hydro minus no-pre mean multiplicity")
    axes[0, 0].legend(frameon=False, fontsize=8)

    for axis, key, label in (
        (axes[0, 1], "delta_zg", r"$\langle z_g^{\rm pre}-z_g^{\rm no}\rangle$"),
        (axes[0, 2], "delta_rg", r"$\langle R_g^{\rm pre}-R_g^{\rm no}\rangle$"),
    ):
        finite_errorbar(
            axis,
            x,
            [all_jets[radius][key].value for radius in RADIUS_DIGITS],
            [all_jets[radius][key].event_error for radius in RADIUS_DIGITS],
            marker="o",
            linestyle="none",
            color="#333333",
            capsize=2.5,
        )
        axis.set_ylabel(label + " (both SD pass)")

    for key, label, marker, color in (
        ("relative_mean_girth_shift", "girth", "o", "#0072B2"),
        ("relative_mean_maxkt_shift", r"maximum $k_T$", "s", "#D55E00"),
    ):
        finite_errorbar(
            axes[1, 0],
            x,
            [100.0 * all_jets[radius][key].value for radius in RADIUS_DIGITS],
            [100.0 * all_jets[radius][key].event_error for radius in RADIUS_DIGITS],
            marker=marker,
            linestyle="none",
            color=color,
            capsize=2.5,
            label=label,
        )
    axes[1, 0].set_ylabel("relative weighted-mean shift [%]")
    axes[1, 0].legend(frameon=False, fontsize=8)

    for key, label, marker, color in (
        ("delta_mean_log10_tauf", "all-split jet mean", "o", "#0072B2"),
        ("delta_hardest_log10_tauf", r"hardest-$k_T$ split", "s", "#D55E00"),
    ):
        finite_errorbar(
            axes[1, 1],
            x,
            [all_jets[radius][key].value for radius in RADIUS_DIGITS],
            [all_jets[radius][key].event_error for radius in RADIUS_DIGITS],
            marker=marker,
            linestyle="none",
            color=color,
            capsize=2.5,
            label=label,
        )
    axes[1, 1].set_ylabel(r"$\Delta\langle\log_{10}(\tau_f/{\rm fm})\rangle$")
    axes[1, 1].legend(frameon=False, fontsize=8)

    finite_errorbar(
        axes[1, 2],
        x,
        [
            100.0 * all_jets[radius]["relative_mean_tauf_shift"].value
            for radius in RADIUS_DIGITS
        ],
        [
            100.0 * all_jets[radius]["relative_mean_tauf_shift"].event_error
            for radius in RADIUS_DIGITS
        ],
        marker="o",
        linestyle="none",
        color="#333333",
        capsize=2.5,
    )
    axes[1, 2].set_ylabel(r"relative arithmetic-mean $\tau_f$ shift [%]")

    for axis in axes.flat:
        axis.axhline(0.0, color="0.55", linewidth=0.9)
        axis.set_xticks(x, labels)
        axis.grid(alpha=0.22)
    figure.suptitle(
        "One-to-one matched all-jet substructure response",
        fontsize=13,
    )
    figure.text(
        0.5,
        0.01,
        "No-prehydro jet selection; matched pre-hydro jet has no pT threshold. "
        "Formation-time entries exist only for R=0.4 and R=0.8.",
        ha="center",
        fontsize=8.5,
    )
    figure.subplots_adjust(
        left=0.075, right=0.985, bottom=0.12, top=0.91, wspace=0.32, hspace=0.3
    )
    stem = out_dir / f"{prefix}_matched_observable_response"
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
        help="maximum pair-axis DeltaR as a fraction of jet radius",
    )
    parser.add_argument(
        "--prefix", default="oo5360_jet_paired_substructure_pt30"
    )
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

    no_suffixes = (
        "Pt",
        "Eta",
        "PtD",
        "NormalPtD",
        "NormalEffectiveMultiplicity",
        "Zg",
        "Rg",
        "Mult",
        "TotalMult",
        "G",
        "MaxKt",
        "SoftDropValid",
        "HardPartonId",
        "PairMatchIndex",
        "PairMatchDR",
        "PairMatchOtherPt",
        "PairMatchOtherHardPartonId",
    )
    pre_suffixes = (
        "Pt",
        "PtD",
        "NormalPtD",
        "Zg",
        "Rg",
        "Mult",
        "TotalMult",
        "G",
        "MaxKt",
        "SoftDropValid",
        "HardPartonId",
    )
    formation_suffixes = (
        "FormationTauF",
        "FormationOffset",
        "FormationHardestValid",
        "FormationHardestTauF",
    )
    with uproot.open(input_root) as root_file:
        pairs = root_file["Pairs"].arrays(
            ["pairId", "eventWeight", "seed", "hydroIndex"], library="np"
        )
        no_branches = ["pairId", "eventWeight"]
        pre_branches = ["pairId", "eventWeight"]
        for radius_digit in RADIUS_DIGITS:
            no_branches.extend(
                f"jet{radius_digit}{suffix}" for suffix in no_suffixes
            )
            pre_branches.extend(
                f"jet{radius_digit}{suffix}" for suffix in pre_suffixes
            )
            if radius_digit in (4, 8):
                no_branches.extend(
                    f"jet{radius_digit}{suffix}" for suffix in formation_suffixes
                )
                pre_branches.extend(
                    f"jet{radius_digit}{suffix}" for suffix in formation_suffixes
                )
        no_arrays = root_file["noPrehydro/Jets"].arrays(
            no_branches, library="ak"
        )
        pre_arrays = root_file["withPrehydro/Jets"].arrays(
            pre_branches, library="ak"
        )

    pair_ids = pairs["pairId"]
    event_weights = pairs["eventWeight"].astype(float)
    hydro_indices = pairs["hydroIndex"].astype(np.int64)
    if len(pair_ids) == 0 or len(np.unique(pair_ids)) != len(pair_ids):
        raise ValueError("Pairs tree has empty or duplicate pairId values")
    if len(np.unique(pairs["seed"])) != len(pair_ids):
        raise ValueError("Pairs tree contains duplicate hard-event seeds")
    if not np.all(np.isfinite(event_weights)) or not np.all(event_weights > 0.0):
        raise ValueError("event weights must be finite and positive")
    for variant, arrays in (
        ("noPrehydro", no_arrays),
        ("withPrehydro", pre_arrays),
    ):
        if not np.array_equal(ak.to_numpy(arrays.pairId), pair_ids):
            raise ValueError(f"{variant} pairId ordering differs from Pairs")
        if not np.array_equal(
            ak.to_numpy(arrays.eventWeight), event_weights
        ):
            raise ValueError(f"{variant} eventWeight differs from Pairs")

    transition_rows: list[list[object]] = []
    fail_rows: list[list[object]] = []
    observable_rows: list[list[object]] = []
    contrast_rows: list[list[object]] = []
    fail_by_radius = {}
    observable_by_radius = {}
    metadata: dict[str, object] = {
        "inputRoot": str(input_root),
        "inputBytes": input_root.stat().st_size,
        "pairCount": len(pair_ids),
        "uniqueSeedCount": len(np.unique(pairs["seed"])),
        "uniqueHydroCount": len(np.unique(hydro_indices)),
        "ptMinGeV": args.pt_min,
        "ptMaxGeV": args.pt_max,
        "absEtaMax": args.abs_eta_max,
        "maxPairMatchDRFraction": args.max_pair_match_dr_fraction,
        "selection": (
            "corrected no-prehydro jet pT/eta and one-to-one pair-axis match; "
            "no matched-jet pT threshold"
        ),
        "softDropTransition": (
            "status in no-prehydro selected jet versus its matched pre-hydro jet"
        ),
        "matchedObservableDefinitions": {
            "multiplicity": "Mult=NNormal+NPositiveWake",
            "totalMultiplicity": "TotalMult=NNormal+NPositiveWake-NNegativeWake",
            "zgRg": "paired differences require SoftDropValid=1 in both jets",
            "formationMeanTauF": (
                "arithmetic mean of all valid exact-angle C/A split tau_f values "
                "within each jet; R=0.4 and R=0.8 only"
            ),
            "formationMeanLog10TauF": (
                "mean log10(tau_f/[fm/c]) over all valid C/A splits within each jet"
            ),
            "formationHardestLog10TauF": (
                "log10(tau_f/[fm/c]) at the global maximum-kT C/A split"
            ),
        },
        "singleCoreCategories": (
            "single-like NormalEffectiveMultiplicity<3; intermediate 3-8; "
            "many-like >=8, all defined on the no-prehydro jet"
        ),
        "weighting": "biased-PYTHIA event weight applied once",
        "eventUncertainty": "paired delete-one-run jackknife",
        "hydroUncertainty": (
            "delete-one-hydroIndex block jackknife retained as a robustness check"
        ),
        "radii": {},
    }

    for radius_digit in RADIUS_DIGITS:
        prefix = f"jet{radius_digit}"
        no_pt = no_arrays[f"{prefix}Pt"]
        no_eta = no_arrays[f"{prefix}Eta"]
        match_index = no_arrays[f"{prefix}PairMatchIndex"]
        match_dr = no_arrays[f"{prefix}PairMatchDR"]
        pre_jet_count = ak.num(pre_arrays[f"{prefix}Pt"], axis=1)
        index_valid = (match_index >= 0) & (match_index < pre_jet_count)
        selection = (
            np.isfinite(no_pt)
            & (no_pt > args.pt_min)
            & np.isfinite(no_eta)
            & (abs(no_eta) < args.abs_eta_max)
            & index_valid
            & np.isfinite(match_dr)
            & (
                match_dr
                < args.max_pair_match_dr_fraction * radius_digit / 10.0
            )
        )
        if args.pt_max is not None:
            selection = selection & (no_pt <= args.pt_max)
        safe_match_index = ak.where(index_valid, match_index, 0)

        for event_matches in ak.to_list(match_index[selection]):
            if len(event_matches) != len(set(event_matches)):
                raise ValueError(
                    f"R=0.{radius_digit} selected matches are not one-to-one"
                )

        selected_lengths = ak.to_numpy(
            ak.sum(selection, axis=1)
        ).astype(np.int64)
        event_indices = np.repeat(
            np.arange(len(event_weights), dtype=np.int64), selected_lengths
        )
        no_values = {
            suffix: selected_flat(
                no_arrays[f"{prefix}{suffix}"], selection
            )
            for suffix in (
                "Pt",
                "PtD",
                "NormalPtD",
                "NormalEffectiveMultiplicity",
                "Zg",
                "Rg",
                "Mult",
                "TotalMult",
                "G",
                "MaxKt",
                "SoftDropValid",
                "HardPartonId",
                "PairMatchOtherPt",
                "PairMatchOtherHardPartonId",
            )
        }
        pre_values = {
            suffix: matched_flat(
                pre_arrays[f"{prefix}{suffix}"],
                safe_match_index,
                selection,
            )
            for suffix in pre_suffixes
        }
        if radius_digit in (4, 8):
            no_mean_tau, no_mean_log_tau = formation_jet_summaries(
                no_arrays[f"{prefix}FormationTauF"],
                no_arrays[f"{prefix}FormationOffset"],
            )
            pre_mean_tau, pre_mean_log_tau = formation_jet_summaries(
                pre_arrays[f"{prefix}FormationTauF"],
                pre_arrays[f"{prefix}FormationOffset"],
            )
            no_values["MeanTauF"] = selected_flat(no_mean_tau, selection)
            no_values["MeanLog10TauF"] = selected_flat(
                no_mean_log_tau, selection
            )
            pre_values["MeanTauF"] = matched_flat(
                pre_mean_tau, safe_match_index, selection
            )
            pre_values["MeanLog10TauF"] = matched_flat(
                pre_mean_log_tau, safe_match_index, selection
            )
            no_values["FormationHardestValid"] = selected_flat(
                no_arrays[f"{prefix}FormationHardestValid"], selection
            )
            no_values["FormationHardestTauF"] = selected_flat(
                no_arrays[f"{prefix}FormationHardestTauF"], selection
            )
            pre_values["FormationHardestValid"] = matched_flat(
                pre_arrays[f"{prefix}FormationHardestValid"],
                safe_match_index,
                selection,
            )
            pre_values["FormationHardestTauF"] = matched_flat(
                pre_arrays[f"{prefix}FormationHardestTauF"],
                safe_match_index,
                selection,
            )
        else:
            empty = np.full(len(event_indices), math.nan)
            no_values["MeanTauF"] = empty.copy()
            no_values["MeanLog10TauF"] = empty.copy()
            pre_values["MeanTauF"] = empty.copy()
            pre_values["MeanLog10TauF"] = empty.copy()
            no_values["FormationHardestValid"] = np.zeros(
                len(event_indices), dtype=int
            )
            no_values["FormationHardestTauF"] = empty.copy()
            pre_values["FormationHardestValid"] = np.zeros(
                len(event_indices), dtype=int
            )
            pre_values["FormationHardestTauF"] = empty.copy()

        no_pt_flat = no_values["Pt"].astype(float)
        pre_pt_flat = pre_values["Pt"].astype(float)
        recorded_other_pt = no_values["PairMatchOtherPt"].astype(float)
        if not np.allclose(
            pre_pt_flat,
            recorded_other_pt,
            rtol=2.0e-6,
            atol=2.0e-5,
            equal_nan=True,
        ):
            raise ValueError(
                f"R=0.{radius_digit} matched pT disagrees with pair audit"
            )
        direct_other_hard_id = pre_values["HardPartonId"].astype(int)
        recorded_other_hard_id = no_values[
            "PairMatchOtherHardPartonId"
        ].astype(int)
        if not np.array_equal(
            direct_other_hard_id, recorded_other_hard_id
        ):
            raise ValueError(
                f"R=0.{radius_digit} hard tag disagrees with pair audit"
            )

        no_sd_status = no_values["SoftDropValid"].astype(int)
        pre_sd_status = pre_values["SoftDropValid"].astype(int)
        if not (
            np.all(np.isin(no_sd_status, (0, 1)))
            and np.all(np.isin(pre_sd_status, (0, 1)))
        ):
            raise ValueError("SoftDropValid must contain only zero or one")
        no_sd_valid = no_sd_status == 1
        pre_sd_valid = pre_sd_status == 1

        epsilon_pre = 1.0 - pre_pt_flat / no_pt_flat
        no_ptd = no_values["PtD"].astype(float)
        pre_ptd = pre_values["PtD"].astype(float)
        no_normal_ptd = no_values["NormalPtD"].astype(float)
        pre_normal_ptd = pre_values["NormalPtD"].astype(float)
        no_zg = no_values["Zg"].astype(float)
        pre_zg = pre_values["Zg"].astype(float)
        no_rg = no_values["Rg"].astype(float)
        pre_rg = pre_values["Rg"].astype(float)
        no_mult = no_values["Mult"].astype(float)
        pre_mult = pre_values["Mult"].astype(float)
        no_total_mult = no_values["TotalMult"].astype(float)
        pre_total_mult = pre_values["TotalMult"].astype(float)
        no_girth = no_values["G"].astype(float)
        pre_girth = pre_values["G"].astype(float)
        no_maxkt = no_values["MaxKt"].astype(float)
        pre_maxkt = pre_values["MaxKt"].astype(float)
        no_mean_tauf = no_values["MeanTauF"].astype(float)
        pre_mean_tauf = pre_values["MeanTauF"].astype(float)
        no_mean_log_tauf = no_values["MeanLog10TauF"].astype(float)
        pre_mean_log_tauf = pre_values["MeanLog10TauF"].astype(float)
        no_hardest_tauf = no_values["FormationHardestTauF"].astype(float)
        pre_hardest_tauf = pre_values["FormationHardestTauF"].astype(float)
        both_sd_valid = no_sd_valid & pre_sd_valid
        both_mean_tauf_valid = (
            np.isfinite(no_mean_tauf)
            & (no_mean_tauf > 0.0)
            & np.isfinite(pre_mean_tauf)
            & (pre_mean_tauf > 0.0)
        )
        both_hardest_tauf_valid = (
            no_values["FormationHardestValid"].astype(int) == 1
        ) & (
            pre_values["FormationHardestValid"].astype(int) == 1
        ) & (
            np.isfinite(no_hardest_tauf)
            & (no_hardest_tauf > 0.0)
            & np.isfinite(pre_hardest_tauf)
            & (pre_hardest_tauf > 0.0)
        )
        observable_values = {
            "epsilon_pre": epsilon_pre,
            "delta_ptd": pre_ptd - no_ptd,
            "relative_ptd": pre_ptd / no_ptd - 1.0,
            "delta_normal_ptd": pre_normal_ptd - no_normal_ptd,
            "relative_normal_ptd": (
                pre_normal_ptd / no_normal_ptd - 1.0
            ),
            "delta_mult": pre_mult - no_mult,
            "delta_total_mult": pre_total_mult - no_total_mult,
            "delta_zg": pre_zg - no_zg,
            "delta_rg": pre_rg - no_rg,
            "delta_girth": pre_girth - no_girth,
            "delta_maxkt": pre_maxkt - no_maxkt,
            "delta_mean_tauf": pre_mean_tauf - no_mean_tauf,
            "delta_mean_log10_tauf": pre_mean_log_tauf - no_mean_log_tauf,
            "delta_hardest_log10_tauf": (
                np.log10(pre_hardest_tauf) - np.log10(no_hardest_tauf)
            ),
        }
        observable_masks = {
            observable: np.ones(len(event_indices), dtype=bool)
            for observable in observable_values
        }
        observable_masks["delta_zg"] = both_sd_valid
        observable_masks["delta_rg"] = both_sd_valid
        observable_masks["delta_mean_tauf"] = both_mean_tauf_valid
        observable_masks["delta_mean_log10_tauf"] = both_mean_tauf_valid
        observable_masks["delta_hardest_log10_tauf"] = both_hardest_tauf_valid
        categories = category_masks(
            no_sd_valid,
            no_values["NormalEffectiveMultiplicity"].astype(float),
        )
        hard_parton_id = no_values["HardPartonId"].astype(int)

        radius_transitions = {}
        radius_fail = {}
        radius_observables = {}
        radius_metadata: dict[str, object] = {
            "selectedMatchedJets": len(event_indices),
            "selectedEvents": int(np.count_nonzero(selected_lengths)),
            "pairPtAudit": "PASS",
            "pairHardTagAudit": "PASS",
            "oneToOneMatchAudit": "PASS",
            "flavors": {},
        }

        for flavor in FLAVORS:
            keep_flavor = flavor_mask(
                hard_parton_id, direct_other_hard_id, flavor
            )
            transition_estimates = {}
            transition_metadata = {}
            for transition, no_pass, pre_pass in TRANSITIONS:
                transition_mask = (
                    keep_flavor
                    & (no_sd_valid == no_pass)
                    & (pre_sd_valid == pre_pass)
                )
                estimate = fraction_estimate(
                    transition_mask,
                    keep_flavor,
                    event_indices,
                    event_weights,
                    hydro_indices,
                )
                audit = mean_estimate(
                    np.ones(len(event_indices), dtype=float),
                    transition_mask,
                    event_indices,
                    event_weights,
                    hydro_indices,
                )
                transition_estimates[transition] = estimate
                transition_metadata[transition] = {
                    **estimate_to_dict(estimate),
                    "numeratorAudit": estimate_to_dict(audit),
                }
                transition_rows.append(
                    [
                        f"0.{radius_digit}",
                        flavor,
                        transition,
                        format_float(estimate.value),
                        format_float(estimate.event_error),
                        format_float(estimate.hydro_error),
                        audit.raw_entries,
                        audit.raw_events,
                        format_float(audit.weighted_entries),
                        estimate.raw_entries,
                        estimate.raw_events,
                        format_float(estimate.weighted_entries),
                        format_float(audit.effective_event_contributors),
                        format_float(audit.max_event_fraction),
                    ]
                )

            no_fail_fraction = fraction_estimate(
                keep_flavor & ~no_sd_valid,
                keep_flavor,
                event_indices,
                event_weights,
                hydro_indices,
            )
            pre_fail_fraction = fraction_estimate(
                keep_flavor & ~pre_sd_valid,
                keep_flavor,
                event_indices,
                event_weights,
                hydro_indices,
            )
            delta_fail_fraction = ratio_estimate(
                numerator=(~pre_sd_valid).astype(float)
                - (~no_sd_valid).astype(float),
                denominator=np.ones(len(event_indices), dtype=float),
                entry_mask=keep_flavor,
                event_indices=event_indices,
                event_weights=event_weights,
                event_hydro_indices=hydro_indices,
            )
            fail_ratio = ratio_estimate(
                numerator=(~pre_sd_valid).astype(float),
                denominator=(~no_sd_valid).astype(float),
                entry_mask=keep_flavor,
                event_indices=event_indices,
                event_weights=event_weights,
                event_hydro_indices=hydro_indices,
            )
            flavor_fail = {
                "no_fail_fraction": no_fail_fraction,
                "pre_fail_fraction": pre_fail_fraction,
                "delta_fail_fraction": delta_fail_fraction,
                "pre_over_no_fail_ratio": fail_ratio,
            }
            for observable, estimate in flavor_fail.items():
                fail_rows.append(
                    [
                        f"0.{radius_digit}",
                        flavor,
                        observable,
                        format_float(estimate.value),
                        format_float(estimate.event_error),
                        format_float(estimate.hydro_error),
                        estimate.raw_entries,
                        estimate.raw_events,
                        format_float(estimate.weighted_entries),
                        format_float(estimate.normalization_denominator),
                        format_float(estimate.effective_event_contributors),
                        format_float(estimate.max_event_fraction),
                    ]
                )

            flavor_observables = {}
            flavor_observable_metadata = {}
            for category, category_mask in categories.items():
                keep_category = keep_flavor & category_mask
                estimates = {
                    observable: mean_estimate(
                        values,
                        keep_category & observable_masks[observable],
                        event_indices,
                        event_weights,
                        hydro_indices,
                    )
                    for observable, values in observable_values.items()
                }
                estimates["pt_weighted_epsilon_pre"] = ratio_estimate(
                    numerator=no_pt_flat * epsilon_pre,
                    denominator=no_pt_flat,
                    entry_mask=keep_category,
                    event_indices=event_indices,
                    event_weights=event_weights,
                    event_hydro_indices=hydro_indices,
                )
                estimates["relative_mean_ptd_shift"] = ratio_estimate(
                    numerator=pre_ptd - no_ptd,
                    denominator=no_ptd,
                    entry_mask=keep_category,
                    event_indices=event_indices,
                    event_weights=event_weights,
                    event_hydro_indices=hydro_indices,
                )
                estimates[
                    "relative_mean_normal_ptd_shift"
                ] = ratio_estimate(
                    numerator=pre_normal_ptd - no_normal_ptd,
                    denominator=no_normal_ptd,
                    entry_mask=keep_category,
                    event_indices=event_indices,
                    event_weights=event_weights,
                    event_hydro_indices=hydro_indices,
                )
                aggregate_relative_inputs = {
                    "relative_mean_mult_shift": (
                        pre_mult - no_mult,
                        no_mult,
                        np.ones(len(event_indices), dtype=bool),
                    ),
                    "relative_mean_total_mult_shift": (
                        pre_total_mult - no_total_mult,
                        no_total_mult,
                        np.ones(len(event_indices), dtype=bool),
                    ),
                    "relative_mean_zg_shift": (
                        pre_zg - no_zg,
                        no_zg,
                        both_sd_valid,
                    ),
                    "relative_mean_rg_shift": (
                        pre_rg - no_rg,
                        no_rg,
                        both_sd_valid,
                    ),
                    "relative_mean_girth_shift": (
                        pre_girth - no_girth,
                        no_girth,
                        np.ones(len(event_indices), dtype=bool),
                    ),
                    "relative_mean_maxkt_shift": (
                        pre_maxkt - no_maxkt,
                        no_maxkt,
                        np.ones(len(event_indices), dtype=bool),
                    ),
                    "relative_mean_tauf_shift": (
                        pre_mean_tauf - no_mean_tauf,
                        no_mean_tauf,
                        both_mean_tauf_valid,
                    ),
                }
                for observable, (
                    numerator,
                    denominator,
                    observable_mask,
                ) in aggregate_relative_inputs.items():
                    estimates[observable] = ratio_estimate(
                        numerator=numerator,
                        denominator=denominator,
                        entry_mask=keep_category & observable_mask,
                        event_indices=event_indices,
                        event_weights=event_weights,
                        event_hydro_indices=hydro_indices,
                    )
                flavor_observables[category] = estimates
                flavor_observable_metadata[category] = {
                    observable: estimate_to_dict(estimate)
                    for observable, estimate in estimates.items()
                }
                for observable in OBSERVABLES:
                    estimate = estimates[observable]
                    observable_rows.append(
                        [
                            f"0.{radius_digit}",
                            flavor,
                            category,
                            observable,
                            format_float(estimate.value),
                            format_float(estimate.event_error),
                            format_float(estimate.hydro_error),
                            estimate.raw_entries,
                            estimate.raw_events,
                            format_float(estimate.weighted_entries),
                            format_float(
                                estimate.normalization_denominator
                            ),
                            format_float(
                                estimate.effective_event_contributors
                            ),
                            format_float(estimate.max_event_fraction),
                        ]
                    )

            single_many_contrasts = {}
            for observable in OBSERVABLES:
                single = flavor_observables["single_like"][observable]
                many = flavor_observables["many_like"][observable]
                contrast = difference_estimate(many, single)
                single_many_contrasts[observable] = difference_to_dict(
                    contrast
                )
                contrast_rows.append(
                    [
                        f"0.{radius_digit}",
                        flavor,
                        observable,
                        format_float(many.value),
                        format_float(single.value),
                        format_float(contrast.value),
                        format_float(contrast.event_error),
                        format_float(contrast.hydro_error),
                        format_float(
                            contrast.value / contrast.event_error
                            if contrast.event_error > 0.0
                            else math.nan
                        ),
                        format_float(
                            contrast.value / contrast.hydro_error
                            if contrast.hydro_error > 0.0
                            else math.nan
                        ),
                        many.raw_entries,
                        single.raw_entries,
                        format_float(many.effective_event_contributors),
                        format_float(single.effective_event_contributors),
                        format_float(many.max_event_fraction),
                        format_float(single.max_event_fraction),
                    ]
                )

            radius_transitions[flavor] = transition_estimates
            radius_fail[flavor] = flavor_fail
            radius_observables[flavor] = flavor_observables
            radius_metadata["flavors"][flavor] = {
                "selectedJets": int(np.count_nonzero(keep_flavor)),
                "selectedEvents": int(
                    np.count_nonzero(
                        np.bincount(
                            event_indices[keep_flavor],
                            minlength=len(event_weights),
                        )
                    )
                ),
                "softDropTransitions": transition_metadata,
                "softDropFailure": {
                    key: estimate_to_dict(value)
                    for key, value in flavor_fail.items()
                },
                "observables": flavor_observable_metadata,
                "manyMinusSingleContrasts": single_many_contrasts,
            }

        transition_sum = sum(
            radius_transitions["all"][transition].value
            for transition, _, _ in TRANSITIONS
        )
        if not math.isfinite(transition_sum) or not math.isclose(
            transition_sum, 1.0, rel_tol=0.0, abs_tol=2.0e-12
        ):
            raise ValueError(
                f"R=0.{radius_digit} Soft Drop transition fractions "
                f"sum to {transition_sum}, not one"
            )

        plot_paths = plot_radius(
            out_dir,
            args.prefix,
            radius_digit,
            args.pt_min,
            args.pt_max,
            radius_transitions,
            radius_fail,
            radius_observables,
        )
        radius_metadata["plots"] = [str(path) for path in plot_paths]
        metadata["radii"][f"0.{radius_digit}"] = radius_metadata
        fail_by_radius[radius_digit] = radius_fail
        observable_by_radius[radius_digit] = radius_observables

    transition_path = out_dir / f"{args.prefix}_softdrop_transitions.tsv"
    with transition_path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "radius",
                "flavor",
                "transition",
                "fraction",
                "eventJackknifeError",
                "hydroBlockJackknifeError",
                "rawTransitionJets",
                "rawTransitionEvents",
                "weightedTransitionJets",
                "rawDenominatorJets",
                "rawDenominatorEvents",
                "weightedDenominatorJets",
                "effectiveTransitionEventContributors",
                "maxTransitionEventContributionFraction",
            ]
        )
        writer.writerows(transition_rows)

    fail_path = out_dir / f"{args.prefix}_softdrop_failure.tsv"
    with fail_path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "radius",
                "flavor",
                "observable",
                "value",
                "eventJackknifeError",
                "hydroBlockJackknifeError",
                "rawEntries",
                "rawEvents",
                "weightedEntries",
                "normalizationDenominator",
                "effectiveEventContributors",
                "maxEventContributionFraction",
            ]
        )
        writer.writerows(fail_rows)

    observable_path = out_dir / f"{args.prefix}_paired_observables.tsv"
    with observable_path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "radius",
                "flavor",
                "category",
                "observable",
                "value",
                "eventJackknifeError",
                "hydroBlockJackknifeError",
                "rawEntries",
                "rawEvents",
                "weightedEntries",
                "normalizationDenominator",
                "effectiveEventContributors",
                "maxEventContributionFraction",
            ]
        )
        writer.writerows(observable_rows)

    contrast_path = out_dir / f"{args.prefix}_single_many_contrasts.tsv"
    with contrast_path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "radius",
                "flavor",
                "observable",
                "manyLikeValue",
                "singleLikeValue",
                "manyMinusSingle",
                "eventJackknifeError",
                "hydroBlockJackknifeError",
                "eventSignificance",
                "hydroSignificance",
                "manyLikeRawEntries",
                "singleLikeRawEntries",
                "manyLikeEffectiveEventContributors",
                "singleLikeEffectiveEventContributors",
                "manyLikeMaxEventContributionFraction",
                "singleLikeMaxEventContributionFraction",
            ]
        )
        writer.writerows(contrast_rows)

    summary_paths = plot_radius_summary(
        out_dir,
        args.prefix,
        fail_by_radius,
        observable_by_radius,
    )
    matched_response_paths = plot_matched_observable_summary(
        out_dir,
        args.prefix,
        observable_by_radius,
    )
    metadata_path = out_dir / f"{args.prefix}_metadata.json"
    metadata["outputFiles"] = {
        "softDropTransitionsTsv": str(transition_path),
        "softDropFailureTsv": str(fail_path),
        "pairedObservablesTsv": str(observable_path),
        "singleManyContrastsTsv": str(contrast_path),
        "summaryPlots": [str(path) for path in summary_paths],
        "matchedObservableResponsePlots": [
            str(path) for path in matched_response_paths
        ],
        "metadataJson": str(metadata_path),
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    return metadata


def main() -> None:
    args = make_parser().parse_args()
    result = run(args)
    print(
        f"Analyzed {result['pairCount']} paired events from "
        f"{result['inputRoot']}"
    )


if __name__ == "__main__":
    main()
