#!/usr/bin/env python3
"""Reconstruct and plot paired OO jet C/A formation-time estimators."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
import numpy as np
import uproot

import plot_oo_jet_variables as base


VARIANTS = ("noPrehydro", "withPrehydro")
VARIANT_LABELS = {
    "noPrehydro": "No pre-hydro",
    "withPrehydro": "Pre-hydro",
}
VARIANT_COLORS = {
    "noPrehydro": "#0072B2",
    "withPrehydro": "#D55E00",
}
VARIANT_MARKERS = {"noPrehydro": "o", "withPrehydro": "s"}
RADII = (0.4, 0.8)
PT_INTERVALS = ((30.0, 50.0), (50.0, 80.0), (80.0, None))
HBARC_GEV_FM = 0.19732698
CURRENT_ROOT_SCHEMA = "oo-paired-root-v8"
ROOT_TAU_F_SCALES = {
    "oo-paired-root-v5": 1.0,
    "oo-paired-root-v6": 1.0,
    "oo-paired-root-v7": 0.5,
    CURRENT_ROOT_SCHEMA: 1.0,
}


def stored_tau_f_scale(schema: str) -> float:
    """Return the scale needed to express stored times in the E/Q^2 convention."""
    try:
        return ROOT_TAU_F_SCALES[schema]
    except KeyError as error:
        raise ValueError(f"unsupported ROOT schema {schema}") from error


def root_string(root_object: Any) -> str:
    """Read a string stored as either a ROOT TNamed or TObjString."""
    if root_object.has_member("fTitle"):
        return str(root_object.member("fTitle"))
    return str(root_object)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def pt_tag(pt_low: float, pt_high: float | None) -> str:
    return f"pt{pt_low:g}to{pt_high:g}" if pt_high is not None else f"pt{pt_low:g}plus"


def pt_label(pt_low: float, pt_high: float | None) -> str:
    if pt_high is None:
        return rf"$p_T^{{\rm jet}}>{pt_low:g}$ GeV"
    return rf"${pt_low:g}<p_T^{{\rm jet}}\leq {pt_high:g}$ GeV"


def jet_selection(
    arrays: ak.Array, radius: float, pt_low: float, pt_high: float | None
) -> np.ndarray:
    selected = (
        np.isclose(ak.to_numpy(arrays.radius), radius, rtol=0.0, atol=1.0e-6)
        & np.isfinite(ak.to_numpy(arrays.correctedPt))
        & (ak.to_numpy(arrays.correctedPt) > pt_low)
    )
    if pt_high is not None:
        selected &= ak.to_numpy(arrays.correctedPt) <= pt_high
    return selected


def flatten_splits(
    arrays: ak.Array, selection: np.ndarray, field: str
) -> tuple[np.ndarray, np.ndarray]:
    selected = arrays[selection]
    values = selected[field]
    lengths = ak.to_numpy(ak.num(values, axis=1)).astype(np.int64)
    flat = ak.to_numpy(ak.flatten(values)).astype(np.float64)
    event_indices = np.repeat(
        ak.to_numpy(selected.eventIndex).astype(np.int64), lengths
    )
    if len(flat) != len(event_indices):
        raise ValueError(f"split/event alignment failed for {field}")
    return flat, event_indices


def flatten_hardest(
    arrays: ak.Array, selection: np.ndarray, field: str
) -> tuple[np.ndarray, np.ndarray]:
    valid = selection & (ak.to_numpy(arrays.hardestValid) == 1)
    values = ak.to_numpy(arrays[field][valid]).astype(np.float64)
    event_indices = ak.to_numpy(arrays.eventIndex[valid]).astype(np.int64)
    finite = np.isfinite(values)
    return values[finite], event_indices[finite]


def event_count(values: np.ndarray, event_count_total: int) -> np.ndarray:
    return np.bincount(values, minlength=event_count_total).astype(np.float64)


def paired_ratio_of_ratios(
    first_numerator: np.ndarray,
    first_denominator: np.ndarray,
    second_numerator: np.ndarray,
    second_denominator: np.ndarray,
    weights: np.ndarray,
) -> tuple[float, float]:
    weighted = [
        weights * counts
        for counts in (
            first_numerator,
            first_denominator,
            second_numerator,
            second_denominator,
        )
    ]
    totals = [float(np.sum(values)) for values in weighted]
    if any(value <= 0.0 for value in totals):
        return math.nan, math.nan
    value = (totals[0] / totals[1]) / (totals[2] / totals[3])
    leave = (
        (totals[0] - weighted[0]) / (totals[1] - weighted[1])
    ) / ((totals[2] - weighted[2]) / (totals[3] - weighted[3]))
    return value, float(base.jackknife_error(leave))


def weighted_histogram_matrix(
    values: np.ndarray,
    event_indices: np.ndarray,
    weights: np.ndarray,
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, int, int, int]:
    finite = np.isfinite(values)
    finite_values = values[finite]
    matrix, underflow, overflow = base.histogram_matrix(
        finite_values, event_indices[finite], weights, edges
    )
    raw_by_bin, _ = np.histogram(finite_values, bins=edges)
    return (
        matrix,
        raw_by_bin.astype(np.int64),
        int(np.count_nonzero(finite)),
        underflow,
        overflow,
    )


def differential_from_matrix(
    matrix: np.ndarray,
    weights: np.ndarray,
    sigma_gen: np.ndarray,
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    widths = np.diff(edges)
    factor, weight_sum, _ = base.pythia_parallel_factor(weights, sigma_gen)
    totals = np.sum(matrix, axis=0)
    values = factor * totals / widths
    weighted_sigma_sum = float(np.sum(weights * sigma_gen))
    leave_values = (
        (weighted_sigma_sum - weights * sigma_gen)[:, np.newaxis]
        * (totals[np.newaxis, :] - matrix)
        / (weight_sum - weights)[:, np.newaxis] ** 2
        / widths[np.newaxis, :]
    )
    return values, base.jackknife_error(leave_values)


def weighted_yield_from_matrix(
    matrix: np.ndarray,
    weights: np.ndarray,
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return event-weighted yield per unit total event weight."""
    widths = np.diff(edges)
    weight_sum = float(np.sum(weights))
    totals = np.sum(matrix, axis=0)
    values = totals / weight_sum / widths
    leave_values = (
        (totals[np.newaxis, :] - matrix)
        / (weight_sum - weights)[:, np.newaxis]
        / widths[np.newaxis, :]
    )
    return values, base.jackknife_error(leave_values)


def selected_jet_normalized_yield(
    matrix: np.ndarray,
    selected_jet_counts: np.ndarray,
    weights: np.ndarray,
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (1/sigma_jet) d sigma_split / dx with paired jackknife errors."""
    widths = np.diff(edges)
    totals = np.sum(matrix, axis=0)
    weighted_jets_by_event = weights * selected_jet_counts
    weighted_jets = float(np.sum(weighted_jets_by_event))
    values = totals / weighted_jets / widths
    leave_values = np.divide(
        totals[np.newaxis, :] - matrix,
        (weighted_jets - weighted_jets_by_event)[:, np.newaxis],
        out=np.full_like(matrix, np.nan),
        where=(weighted_jets - weighted_jets_by_event)[:, np.newaxis] != 0.0,
    ) / widths[np.newaxis, :]
    return values, base.jackknife_error(leave_values)


def weighted_quantile(
    values: np.ndarray, sample_weights: np.ndarray, quantiles: list[float]
) -> list[float]:
    finite = np.isfinite(values) & np.isfinite(sample_weights) & (sample_weights > 0.0)
    values = values[finite]
    sample_weights = sample_weights[finite]
    if len(values) == 0:
        return [math.nan for _ in quantiles]
    order = np.argsort(values)
    values = values[order]
    sample_weights = sample_weights[order]
    cumulative = np.cumsum(sample_weights) - 0.5 * sample_weights
    cumulative /= np.sum(sample_weights)
    return [float(np.interp(value, cumulative, values)) for value in quantiles]


def finite_errorbar(axis, x, y, yerr, **kwargs) -> None:
    x = np.asarray(x)
    y = np.asarray(y)
    yerr = np.asarray(yerr)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(yerr)
    if np.any(valid):
        axis.errorbar(x[valid], y[valid], yerr=yerr[valid], **kwargs)


def spectrum_plot(
    results: dict[str, dict[str, dict[str, Any]]],
    radius: float,
    pt_low: float,
    pt_high: float | None,
    edges: np.ndarray,
    output_base: Path,
    sample_label: str,
) -> None:
    centers = 0.5 * (edges[:-1] + edges[1:])
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(11.8, 7.6),
        sharex="col",
        gridspec_kw={"height_ratios": [3.0, 1.15]},
    )
    for column, split_kind in enumerate(("all", "hardest_kt")):
        for variant in VARIANTS:
            item = results[split_kind][variant]
            finite_errorbar(
                axes[0, column],
                centers,
                item["normalized_values"],
                item["normalized_errors"],
                color=VARIANT_COLORS[variant],
                marker=VARIANT_MARKERS[variant],
                markersize=3.2,
                linewidth=1.0,
                capsize=1.8,
                label=VARIANT_LABELS[variant],
            )
        ratio = results[split_kind]["ratio"]
        finite_errorbar(
            axes[1, column],
            centers,
            ratio["plot_values"],
            ratio["plot_errors"],
            color="#333333",
            marker="o",
            markersize=3.2,
            linewidth=1.0,
            capsize=1.8,
        )
        axes[0, column].set_yscale("log")
        axes[0, column].set_ylabel(
            r"$(1/\sigma_{\rm jet})\,d\sigma_{\rm split}/d\log_{10}\tau_{\rm f}$"
        )
        axes[1, column].set_ylabel("Pre-hydro / no pre-hydro")
        axes[1, column].set_xlabel(r"$\log_{10}(\tau_{\rm f}/[\mathrm{fm}/c])$")
        axes[1, column].axhline(1.0, color="0.5", linewidth=0.9)
        axes[1, column].set_ylim(
            base.ratio_limits(ratio["plot_values"], ratio["plot_errors"])
        )
        for row in range(2):
            axes[row, column].grid(alpha=0.2)
            axes[row, column].set_xlim(edges[0], edges[-1])
            axes[row, column].axvline(math.log10(0.1), color="#009E73", linestyle="--", linewidth=0.9)
            axes[row, column].axvline(math.log10(0.24), color="#CC79A7", linestyle=":", linewidth=1.0)
        axes[0, column].set_title(
            "All C/A declusterings"
            if split_kind == "all"
            else r"global hardest-$k_T$ split per jet"
        )
    axes[0, 0].legend(frameon=False, fontsize=9)
    figure.suptitle(
        rf"O+O 5.36 TeV, 0--5%, anti-$k_T$ R={radius:g}, {pt_label(pt_low, pt_high)}"
        + f"\n{sample_label}; positive-constituent C/A tree; "
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


def correlation_axes(radius: float) -> dict[str, tuple[np.ndarray, str, str]]:
    return {
        "z": (np.linspace(0.0, 0.5, 21), "z", r"$z=\min(E_1,E_2)/E_{\rm parent}$"),
        "deltaR": (
            np.linspace(-3.0, math.log10(2.0 * radius + 1.0e-6), 21),
            "log_delta_r",
            r"$\log_{10}\Delta R_{12}$",
        ),
        "kt": (
            np.linspace(-3.0, 3.0, 25),
            "log_kt",
            r"$\log_{10}(k_T/\mathrm{GeV})$",
        ),
    }


def transform_correlation(name: str, values: np.ndarray) -> np.ndarray:
    if name == "z":
        return values
    transformed = np.full_like(values, np.nan, dtype=np.float64)
    positive = values > 0.0
    transformed[positive] = np.log10(values[positive])
    return transformed


def correlation_plot(
    correlation_data: dict[str, dict[str, Any]],
    radius: float,
    pt_low: float,
    pt_high: float | None,
    tau_edges: np.ndarray,
    output_base: Path,
    sample_label: str,
) -> None:
    specs = correlation_axes(radius)
    figure, axes = plt.subplots(3, 3, figsize=(13.6, 10.0), sharey=True)
    column_titles = ("No pre-hydro", "Pre-hydro", "Pre-hydro / no pre-hydro")
    for column, title in enumerate(column_titles):
        axes[0, column].set_title(title)
    for row, (name, (_, _, x_label)) in enumerate(specs.items()):
        no_density = correlation_data[name]["noPrehydro"]["density"]
        pre_density = correlation_data[name]["withPrehydro"]["density"]
        x_edges = correlation_data[name]["x_edges"]
        positive = np.concatenate((no_density[no_density > 0.0], pre_density[pre_density > 0.0]))
        if len(positive) == 0:
            vmin, vmax = 1.0e-12, 1.0
        else:
            vmin = max(float(np.percentile(positive, 3.0)), float(np.min(positive)))
            vmax = float(np.max(positive))
            if vmin >= vmax:
                vmin = vmax / 10.0
        for column, variant in enumerate(VARIANTS):
            density = np.ma.masked_less_equal(
                correlation_data[name][variant]["density"].T, 0.0
            )
            mesh = axes[row, column].pcolormesh(
                x_edges,
                tau_edges,
                density,
                shading="auto",
                cmap="viridis",
                norm=LogNorm(vmin=vmin, vmax=vmax),
            )
            figure.colorbar(
                mesh,
                ax=axes[row, column],
                pad=0.01,
                label=r"$(1/\sigma_{\rm jet})d^2\sigma_{\rm split}/dx\,d\log\tau_f$",
            )
        ratio = np.ma.masked_invalid(correlation_data[name]["ratio"].T)
        ratio = np.ma.masked_where(correlation_data[name]["no_raw"].T < 10, ratio)
        finite_ratio = ratio.compressed()
        if len(finite_ratio):
            distance = max(
                0.2,
                min(1.0, float(np.nanpercentile(np.abs(finite_ratio - 1.0), 95.0))),
            )
        else:
            distance = 0.5
        mesh = axes[row, 2].pcolormesh(
            x_edges,
            tau_edges,
            ratio,
            shading="auto",
            cmap="coolwarm",
            norm=TwoSlopeNorm(vmin=1.0 - distance, vcenter=1.0, vmax=1.0 + distance),
        )
        figure.colorbar(mesh, ax=axes[row, 2], pad=0.01, label="ratio")
        for column in range(3):
            axes[row, column].set_xlabel(x_label)
            axes[row, column].set_ylabel(r"$\log_{10}(\tau_{\rm f}/[\mathrm{fm}/c])$")
    figure.suptitle(
        rf"All C/A declusterings: anti-$k_T$ R={radius:g}, {pt_label(pt_low, pt_high)}"
        + f"\n{sample_label}; ratio cells require at least 10 no-prehydro entries",
        fontsize=11,
    )
    figure.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    for suffix in ("pdf", "png"):
        figure.savefig(output_base.with_suffix(f".{suffix}"), dpi=180)
    plt.close(figure)


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def radius_piece(
    events: dict[str, np.ndarray],
    arrays: ak.Array,
    radius_digit: int,
    tau_f_scale: float,
) -> ak.Array:
    prefix = f"jet{radius_digit}"
    pt = arrays[f"{prefix}Pt"]
    jet_counts = ak.to_numpy(ak.num(pt, axis=1)).astype(np.int64)
    event_indices = np.repeat(
        np.arange(len(events["pairId"]), dtype=np.int64), jet_counts
    )

    def flat(suffix: str) -> ak.Array:
        return ak.flatten(arrays[f"{prefix}{suffix}"], axis=1)

    offsets = arrays[f"{prefix}FormationOffset"]
    if not bool(
        ak.all(ak.num(offsets, axis=1) == ak.num(pt, axis=1) + 1)
        and ak.all(ak.firsts(offsets) == 0)
        and ak.all(
            offsets[:, -1]
            == ak.num(arrays[f"{prefix}FormationTauF"], axis=1)
        )
    ):
        raise ValueError(f"{prefix} formation-time offsets are inconsistent")
    split_counts = ak.to_numpy(
        ak.flatten(offsets[:, 1:] - offsets[:, :-1], axis=1)
    ).astype(np.int64)

    def split_flat(suffix: str) -> ak.Array:
        values = ak.flatten(arrays[f"{prefix}{suffix}"], axis=1)
        return ak.unflatten(values, split_counts)

    return ak.zip(
        {
            "pairId": events["pairId"][event_indices],
            "eventIndex": event_indices,
            "eventWeight": events["eventWeight"][event_indices],
            "sigmaGen": events["sigmaGen"][event_indices],
            "radius": np.full(len(event_indices), radius_digit / 10.0),
            "jetIndex": ak.to_numpy(ak.flatten(ak.local_index(pt, axis=1))),
            "correctedPt": flat("Pt"),
            "correctedEta": flat("Eta"),
            "correctedPhi": flat("Phi"),
            "rawPt": flat("RawPt"),
            "multiplicity": flat("Mult"),
            "nNormal": flat("NNormal"),
            "nPositiveWake": flat("NPositiveWake"),
            "nNegativeWake": flat("NNegativeWake"),
            "nNegativeThermal": flat("NNegativeThermal"),
            "nHadronizedHoles": flat("NHadronizedHoles"),
            "positiveWakePt": flat("PositiveWakePt"),
            "negativeWakePt": flat("NegativeWakePt"),
            "nInvalidSplittings": flat("FormationInvalidSplits"),
            "tauF": tau_f_scale * split_flat("FormationTauF"),
            "tauFSmallAngle": tau_f_scale * split_flat("FormationTauFSmallAngle"),
            "z": split_flat("FormationZ"),
            "theta": split_flat("FormationTheta"),
            "deltaR": split_flat("FormationDeltaR"),
            "kt": split_flat("FormationKt"),
            "parentEnergy": split_flat("FormationParentE"),
            "hardestValid": flat("FormationHardestValid"),
            "hardestTauF": tau_f_scale * flat("FormationHardestTauF"),
            "hardestTauFSmallAngle": (
                tau_f_scale * flat("FormationHardestTauFSmallAngle")
            ),
            "hardestZ": flat("FormationHardestZ"),
            "hardestTheta": flat("FormationHardestTheta"),
            "hardestDeltaR": flat("FormationHardestDeltaR"),
            "hardestKt": flat("FormationHardestKt"),
            "hardestParentEnergy": flat("FormationHardestParentE"),
        },
        depth_limit=1,
    )


def validate_root(
    input_root: Path,
) -> tuple[dict[str, np.ndarray], dict[str, ak.Array], str, float]:
    jet_suffixes = (
        "Pt",
        "Eta",
        "Phi",
        "RawPt",
        "Mult",
        "NNormal",
        "NPositiveWake",
        "NNegativeWake",
        "NNegativeThermal",
        "NHadronizedHoles",
        "PositiveWakePt",
        "NegativeWakePt",
        "FormationInvalidSplits",
        "FormationTauF",
        "FormationTauFSmallAngle",
        "FormationZ",
        "FormationTheta",
        "FormationDeltaR",
        "FormationKt",
        "FormationParentE",
        "FormationOffset",
        "FormationHardestValid",
        "FormationHardestTauF",
        "FormationHardestTauFSmallAngle",
        "FormationHardestZ",
        "FormationHardestTheta",
        "FormationHardestDeltaR",
        "FormationHardestKt",
        "FormationHardestParentE",
    )
    branches = ["pairId", "eventWeight", "sigmaGen"]
    for radius_digit in (4, 8):
        branches.extend(f"jet{radius_digit}{suffix}" for suffix in jet_suffixes)
    with uproot.open(input_root) as root_file:
        schema = root_string(root_file["metadata/schemaVersion"])
        if schema not in ROOT_TAU_F_SCALES:
            raise ValueError(
                "formation-time analysis requires oo-paired-root-v5 through v8; "
                f"found {schema}"
            )
        tau_f_scale = stored_tau_f_scale(schema)
        events = root_file["Pairs"].arrays(
            ["pairId", "eventWeight", "sigmaGen", "seed", "hydroIndex"],
            library="np",
        )
        variants: dict[str, ak.Array] = {}
        for variant in VARIANTS:
            arrays = root_file[f"{variant}/Jets"].arrays(branches, library="ak")
            if not np.array_equal(ak.to_numpy(arrays.pairId), events["pairId"]):
                raise ValueError(f"{variant} pair IDs disagree with Pairs")
            if not np.array_equal(
                ak.to_numpy(arrays.eventWeight), events["eventWeight"]
            ) or not np.array_equal(
                ak.to_numpy(arrays.sigmaGen), events["sigmaGen"]
            ):
                raise ValueError(f"{variant} event normalization disagrees with Pairs")
            variants[variant] = ak.concatenate(
                [
                    radius_piece(events, arrays, 4, tau_f_scale),
                    radius_piece(events, arrays, 8, tau_f_scale),
                ],
                axis=0,
            )
    event_count_total = len(events["pairId"])
    if event_count_total == 0 or len(np.unique(events["pairId"])) != event_count_total:
        raise ValueError("Events tree is empty or has duplicate pair IDs")
    weights = events["eventWeight"].astype(np.float64)
    sigma_gen = events["sigmaGen"].astype(np.float64)
    if not np.all(np.isfinite(weights) & (weights > 0.0)):
        raise ValueError("event weights must be finite and positive")
    if not np.all(np.isfinite(sigma_gen) & (sigma_gen > 0.0)):
        raise ValueError("sigmaGen values must be finite and positive")
    for variant, arrays in variants.items():
        indices = ak.to_numpy(arrays.eventIndex).astype(np.int64)
        if np.any((indices < 0) | (indices >= event_count_total)):
            raise ValueError(f"{variant} has an invalid eventIndex")
        if not np.array_equal(
            ak.to_numpy(arrays.pairId), events["pairId"][indices]
        ):
            raise ValueError(f"{variant} pair IDs disagree with Events")
        if not np.array_equal(
            ak.to_numpy(arrays.eventWeight), weights[indices]
        ) or not np.array_equal(ak.to_numpy(arrays.sigmaGen), sigma_gen[indices]):
            raise ValueError(f"{variant} event normalization disagrees with Events")
        radius_values = ak.to_numpy(arrays.radius).astype(np.float64)
        supported_radius = np.any(
            np.isclose(radius_values[:, np.newaxis], np.asarray(RADII), atol=1.0e-6),
            axis=1,
        )
        if not np.all(supported_radius):
            raise ValueError(f"{variant} contains an unsupported radius")
        split_lengths = ak.num(arrays.tauF, axis=1)
        for field in (
            "tauFSmallAngle",
            "z",
            "theta",
            "deltaR",
            "kt",
            "parentEnergy",
        ):
            if not bool(ak.all(ak.num(arrays[field], axis=1) == split_lengths)):
                raise ValueError(f"{variant} split-vector length mismatch for {field}")
        retained = np.isfinite(ak.to_numpy(arrays.correctedPt)) & (
            ak.to_numpy(arrays.correctedPt) > 30.0
        )
        valid_count = ak.to_numpy(split_lengths).astype(np.int64)
        invalid_count = ak.to_numpy(arrays.nInvalidSplittings).astype(np.int64)
        multiplicity = ak.to_numpy(arrays.multiplicity).astype(np.int64)
        expected_count = np.maximum(multiplicity - 1, 0)
        if np.any(
            valid_count[retained] + invalid_count[retained]
            != expected_count[retained]
        ):
            raise ValueError(f"{variant} C/A tree does not close at Nconstituent-1")
        if np.any(valid_count[~retained] != 0) or np.any(
            invalid_count[~retained] != 0
        ):
            raise ValueError(
                f"{variant} stores formation splits outside the declared scope"
            )
        hardest_valid = ak.to_numpy(arrays.hardestValid).astype(np.int64) == 1
        if not np.array_equal(hardest_valid, valid_count > 0):
            raise ValueError(f"{variant} hardest-split validity is inconsistent")
        flat_tau = ak.to_numpy(ak.flatten(arrays.tauF)).astype(np.float64)
        flat_tau_small = ak.to_numpy(ak.flatten(arrays.tauFSmallAngle)).astype(
            np.float64
        )
        flat_z = ak.to_numpy(ak.flatten(arrays.z)).astype(np.float64)
        flat_theta = ak.to_numpy(ak.flatten(arrays.theta)).astype(np.float64)
        flat_delta_r = ak.to_numpy(ak.flatten(arrays.deltaR)).astype(np.float64)
        flat_kt = ak.to_numpy(ak.flatten(arrays.kt)).astype(np.float64)
        flat_parent_energy = ak.to_numpy(ak.flatten(arrays.parentEnergy)).astype(
            np.float64
        )
        physical = (
            np.isfinite(flat_tau)
            & (flat_tau > 0.0)
            & np.isfinite(flat_tau_small)
            & (flat_tau_small > 0.0)
            & np.isfinite(flat_z)
            & (flat_z > 0.0)
            & (flat_z <= 0.5 + 1.0e-6)
            & np.isfinite(flat_theta)
            & (flat_theta > 0.0)
            & (flat_theta <= math.pi)
            & np.isfinite(flat_delta_r)
            & (flat_delta_r > 0.0)
            & np.isfinite(flat_kt)
            & (flat_kt > 0.0)
            & np.isfinite(flat_parent_energy)
            & (flat_parent_energy > 0.0)
        )
        if not np.all(physical):
            raise ValueError(f"{variant} has a nonphysical stored formation split")
        valid_hardest = arrays.hardestValid == 1
        if not bool(
            ak.all(
                abs(ak.max(arrays.kt, axis=1, mask_identity=True)[valid_hardest]
                    - arrays.hardestKt[valid_hardest])
                < 2.0e-4
            )
        ):
            raise ValueError(f"{variant} hardest-kT audit failed")
    return events, variants, schema, tau_f_scale


def analyze(args: argparse.Namespace) -> dict[str, Any]:
    input_root = args.input_root.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    events, variants, input_root_schema, input_tau_f_scale = validate_root(
        input_root
    )
    weights = events["eventWeight"].astype(np.float64)
    sigma_gen = events["sigmaGen"].astype(np.float64)
    event_count_total = len(weights)
    factor, weight_sum, sigma_merged = base.pythia_parallel_factor(
        weights, sigma_gen
    )
    effective_events = float(np.sum(weights) ** 2 / np.sum(weights**2))
    tau_edges = np.linspace(args.log_tau_min, args.log_tau_max, args.log_tau_bins + 1)

    spectra_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    correlation_rows: list[dict[str, Any]] = []
    metadata: dict[str, Any] = {
        "schemaVersion": "oo-jet-formation-time-v3",
        "inputRoot": str(input_root),
        "inputRootSchema": input_root_schema,
        "storedTauFScaleApplied": input_tau_f_scale,
        "inputRootBytes": input_root.stat().st_size,
        "inputRootSha256": sha256(input_root),
        "pairCount": event_count_total,
        "weightSum": weight_sum,
        "sigmaMergedMb": sigma_merged,
        "pythiaParallelFactorMb": factor,
        "effectiveEventCount": effective_events,
        "radii": list(RADII),
        "correctedPtIntervalsGeV": [
            [low, high] for low, high in PT_INTERVALS
        ],
        "log10TauEdges": tau_edges.tolist(),
        "hbarCGeVFm": HBARC_GEV_FM,
        "formationTime": (
            "hbarc*Eparent/Qparent^2 = "
            "hbarc/[2*Eparent*z1*z2*(1-cos(theta12))], zi=Ei/Eparent; "
            "massless daughters and exact three-dimensional opening angle"
        ),
        "coefficientConvention": (
            "E/Q^2 convention; schema-v7 2E/Q^2 values are divided by two on read"
        ),
        "smallAngleValidation": (
            "hbarc/[Eparent*z1*z2*DeltaR12^2]; retained only as a validation"
        ),
        "allSplittings": "every valid internal node in the full recursive C/A tree",
        "singleSplitting": (
            "global maximum of min(pT1,pT2)*DeltaR12 over the full C/A tree; "
            "one entry per selected jet when a valid split exists"
        ),
        "zForCorrelations": "min(E1,E2)/Eparent",
        "ktDefinition": "min(pT1,pT2)*DeltaR12",
        "jetSelection": (
            "independent in each variant using 4MomSub-corrected jet pT; lower edge "
            "exclusive, upper edge inclusive; no additional eta cut beyond ROOT storage"
        ),
        "constituentTree": (
            "normal raw-label-0 plus positive-wake raw-label-1 hadrons at physical "
            "four-momentum, C/A E-scheme with R_CA=2*R_antiKt+1e-6 and "
            "FastJet Best strategy"
        ),
        "negativeTreatment": (
            "negative-wake raw-label-2 particles and raw-label-3 hadronized holes are "
            "ghost-associated and subtracted only at jet four-vector level; they set the "
            "corrected-pT selection but are excluded from the nonlinear C/A tree"
        ),
        "normalization": (
            "plotted (1/sigmaJetSelected)*dSigmaSplit/dlog10(tau) is the weighted "
            "split count divided by the weighted selected-jet count and bin width; "
            "event-normalized yield and absolute dSigma/dlog10(tau) remain tabulated; "
            "the selected-jet denominator is recomputed in every jackknife replica"
        ),
        "uncertainty": "paired delete-one-run jackknife",
        "ratioPlotMinimumNoPrehydroEntriesPerBin": 20,
        "interpretation": (
            "Cambridge-Aachen declustering is a final-constituent formation-time "
            "estimator and does not reproduce generator-level parton-shower history"
        ),
        "selections": {},
    }

    for radius in RADII:
        for pt_low, pt_high in PT_INTERVALS:
            tag = f"R{int(round(10 * radius)):02d}_{pt_tag(pt_low, pt_high)}"
            selections = {
                variant: jet_selection(variants[variant], radius, pt_low, pt_high)
                for variant in VARIANTS
            }
            result: dict[str, dict[str, dict[str, Any]]] = {
                "all": {},
                "hardest_kt": {},
            }
            selection_metadata: dict[str, Any] = {}
            integrated_counts: dict[str, dict[str, np.ndarray]] = {}
            selection_summary_rows: dict[str, dict[str, Any]] = {}
            for variant in VARIANTS:
                arrays = variants[variant]
                selection = selections[variant]
                selected_indices = ak.to_numpy(arrays.eventIndex[selection]).astype(
                    np.int64
                )
                selected_jet_counts = event_count(selected_indices, event_count_total)
                integrated_jets = base.integrated_cross_section(
                    selected_jet_counts, weights, sigma_gen
                )
                invalid_splits = int(ak.sum(arrays.nInvalidSplittings[selection]))
                all_tau, all_event_indices = flatten_splits(
                    arrays, selection, "tauF"
                )
                all_tau_small, small_event_indices = flatten_splits(
                    arrays, selection, "tauFSmallAngle"
                )
                if not np.array_equal(all_event_indices, small_event_indices):
                    raise ValueError("exact/small-angle split ordering differs")
                all_log_tau = np.log10(all_tau)
                (
                    all_matrix,
                    raw_all_by_bin,
                    raw_all,
                    all_underflow,
                    all_overflow,
                ) = weighted_histogram_matrix(
                    all_log_tau, all_event_indices, weights, tau_edges
                )
                all_values, all_errors = weighted_yield_from_matrix(
                    all_matrix, weights, tau_edges
                )
                all_normalized, all_normalized_errors = selected_jet_normalized_yield(
                    all_matrix, selected_jet_counts, weights, tau_edges
                )
                all_cross_section, all_cross_section_errors = differential_from_matrix(
                    all_matrix, weights, sigma_gen, tau_edges
                )
                hardest_tau, hardest_event_indices = flatten_hardest(
                    arrays, selection, "hardestTauF"
                )
                integrated_counts[variant] = {
                    "selected_jets": selected_jet_counts,
                    "all_splits": event_count(all_event_indices, event_count_total),
                    "hardest_splits": event_count(
                        hardest_event_indices, event_count_total
                    ),
                }
                hardest_log_tau = np.log10(hardest_tau)
                (
                    hardest_matrix,
                    raw_hardest_by_bin,
                    raw_hardest,
                    hardest_underflow,
                    hardest_overflow,
                ) = weighted_histogram_matrix(
                    hardest_log_tau, hardest_event_indices, weights, tau_edges
                )
                hardest_values, hardest_errors = weighted_yield_from_matrix(
                    hardest_matrix, weights, tau_edges
                )
                hardest_normalized, hardest_normalized_errors = (
                    selected_jet_normalized_yield(
                        hardest_matrix, selected_jet_counts, weights, tau_edges
                    )
                )
                hardest_cross_section, hardest_cross_section_errors = differential_from_matrix(
                    hardest_matrix, weights, sigma_gen, tau_edges
                )
                result["all"][variant] = {
                    "matrix": all_matrix,
                    "values": all_values,
                    "errors": all_errors,
                    "normalized_values": all_normalized,
                    "normalized_errors": all_normalized_errors,
                    "cross_section": all_cross_section,
                    "cross_section_errors": all_cross_section_errors,
                    "raw": raw_all,
                    "raw_by_bin": raw_all_by_bin,
                    "underflow": all_underflow,
                    "overflow": all_overflow,
                }
                result["hardest_kt"][variant] = {
                    "matrix": hardest_matrix,
                    "values": hardest_values,
                    "errors": hardest_errors,
                    "normalized_values": hardest_normalized,
                    "normalized_errors": hardest_normalized_errors,
                    "cross_section": hardest_cross_section,
                    "cross_section_errors": hardest_cross_section_errors,
                    "raw": raw_hardest,
                    "raw_by_bin": raw_hardest_by_bin,
                    "underflow": hardest_underflow,
                    "overflow": hardest_overflow,
                }

                split_ratio = all_tau_small / all_tau
                split_weights = weights[all_event_indices]
                p16_ratio, median_ratio, p84_ratio = weighted_quantile(
                    split_ratio,
                    split_weights,
                    [0.16, 0.5, 0.84],
                )
                p95_abs_log = weighted_quantile(
                    np.abs(np.log10(split_ratio)), split_weights, [0.95]
                )[0]
                weighted_jets = float(np.sum(weights * selected_jet_counts))
                weighted_splits = float(np.sum(weights[all_event_indices]))
                weighted_hardest = float(np.sum(weights[hardest_event_indices]))
                mean_splits_per_jet = (
                    weighted_splits / weighted_jets if weighted_jets > 0.0 else math.nan
                )
                positive_wake_count = event_count(
                    ak.to_numpy(arrays.eventIndex[selection & (ak.to_numpy(arrays.nPositiveWake) > 0)]).astype(np.int64),
                    event_count_total,
                )
                negative_wake_count = event_count(
                    ak.to_numpy(arrays.eventIndex[selection & (ak.to_numpy(arrays.nNegativeWake) > 0)]).astype(np.int64),
                    event_count_total,
                )
                positive_fraction = base.paired_integrated_ratio(
                    positive_wake_count, selected_jet_counts, weights
                )[0] if weighted_jets > 0.0 else math.nan
                negative_fraction = base.paired_integrated_ratio(
                    negative_wake_count, selected_jet_counts, weights
                )[0] if weighted_jets > 0.0 else math.nan
                summary_row = {
                    "variant": variant,
                    "radius": f"{radius:.1f}",
                    "pt_low_exclusive_GeV": f"{pt_low:g}",
                    "pt_high_inclusive_GeV": (
                        "inf" if pt_high is None else f"{pt_high:g}"
                    ),
                    "raw_selected_jets": int(np.sum(selected_jet_counts)),
                    "events_with_selected_jets": int(
                        np.count_nonzero(selected_jet_counts)
                    ),
                    "weighted_selected_jets": f"{weighted_jets:.12e}",
                    "selected_jet_cross_section_mb": (
                        f"{integrated_jets.cross_section:.12e}"
                    ),
                    "selected_jet_cross_section_error_mb": (
                        f"{integrated_jets.error:.12e}"
                    ),
                    "raw_valid_all_splits": raw_all,
                    "raw_valid_hardest_splits": raw_hardest,
                    "invalid_splits": invalid_splits,
                    "weighted_all_splits": f"{weighted_splits:.12e}",
                    "weighted_hardest_splits": f"{weighted_hardest:.12e}",
                    "weighted_splits_per_selected_jet": (
                        f"{mean_splits_per_jet:.12e}"
                    ),
                    "weighted_positive_wake_jet_fraction": (
                        f"{positive_fraction:.12e}"
                    ),
                    "weighted_negative_wake_or_hole_jet_fraction": (
                        f"{negative_fraction:.12e}"
                    ),
                    "small_over_exact_tau_median": f"{median_ratio:.12e}",
                    "small_over_exact_tau_p16": f"{p16_ratio:.12e}",
                    "small_over_exact_tau_p84": f"{p84_ratio:.12e}",
                    "abs_log10_small_over_exact_tau_p95": (
                        f"{p95_abs_log:.12e}"
                    ),
                }
                summary_rows.append(summary_row)
                selection_summary_rows[variant] = summary_row
                selection_metadata[variant] = {
                    "rawSelectedJets": int(np.sum(selected_jet_counts)),
                    "eventsWithSelectedJets": int(np.count_nonzero(selected_jet_counts)),
                    "rawValidAllSplits": raw_all,
                    "rawValidHardestSplits": raw_hardest,
                    "invalidSplits": invalid_splits,
                    "weightedSplitsPerSelectedJet": mean_splits_per_jet,
                    "weightedPositiveWakeJetFraction": positive_fraction,
                    "weightedNegativeWakeOrHoleJetFraction": negative_fraction,
                    "smallOverExactTauMedian": median_ratio,
                    "smallOverExactTauP16": p16_ratio,
                    "smallOverExactTauP84": p84_ratio,
                    "absLog10SmallOverExactTauP95": p95_abs_log,
                }

            for count_name, field_name in (
                ("selected_jets", "selected_jet_yield"),
                ("all_splits", "all_split_yield"),
                ("hardest_splits", "hardest_split_yield"),
            ):
                integrated_ratio, integrated_error = base.paired_integrated_ratio(
                    integrated_counts["withPrehydro"][count_name],
                    integrated_counts["noPrehydro"][count_name],
                    weights,
                )
                selection_summary_rows["noPrehydro"][
                    f"pre_over_no_{field_name}"
                ] = "1.000000000000e+00"
                selection_summary_rows["noPrehydro"][
                    f"pre_over_no_{field_name}_stat_error"
                ] = "0.000000000000e+00"
                selection_summary_rows["withPrehydro"][
                    f"pre_over_no_{field_name}"
                ] = f"{integrated_ratio:.12e}"
                selection_summary_rows["withPrehydro"][
                    f"pre_over_no_{field_name}_stat_error"
                ] = f"{integrated_error:.12e}"
                selection_metadata["withPrehydro"][
                    f"preOverNo{field_name.title().replace('_', '')}"
                ] = integrated_ratio
                selection_metadata["withPrehydro"][
                    f"preOverNo{field_name.title().replace('_', '')}StatError"
                ] = integrated_error
            split_rate_ratio, split_rate_error = paired_ratio_of_ratios(
                integrated_counts["withPrehydro"]["all_splits"],
                integrated_counts["withPrehydro"]["selected_jets"],
                integrated_counts["noPrehydro"]["all_splits"],
                integrated_counts["noPrehydro"]["selected_jets"],
                weights,
            )
            for variant, value, error in (
                ("noPrehydro", 1.0, 0.0),
                ("withPrehydro", split_rate_ratio, split_rate_error),
            ):
                selection_summary_rows[variant][
                    "pre_over_no_splits_per_selected_jet"
                ] = f"{value:.12e}"
                selection_summary_rows[variant][
                    "pre_over_no_splits_per_selected_jet_stat_error"
                ] = f"{error:.12e}"
            selection_metadata["withPrehydro"][
                "preOverNoSplitsPerSelectedJet"
            ] = split_rate_ratio
            selection_metadata["withPrehydro"][
                "preOverNoSplitsPerSelectedJetStatError"
            ] = split_rate_error

            for split_kind in ("all", "hardest_kt"):
                absolute_ratio, absolute_ratio_error = base.paired_ratio(
                    result[split_kind]["withPrehydro"]["matrix"],
                    result[split_kind]["noPrehydro"]["matrix"],
                )
                ratio, ratio_error = base.paired_normalized_ratio(
                    result[split_kind]["withPrehydro"]["matrix"],
                    result[split_kind]["noPrehydro"]["matrix"],
                    integrated_counts["withPrehydro"]["selected_jets"],
                    integrated_counts["noPrehydro"]["selected_jets"],
                    weights,
                )
                result[split_kind]["ratio"] = {
                    "values": ratio,
                    "errors": ratio_error,
                    "plot_values": np.where(
                        result[split_kind]["noPrehydro"]["raw_by_bin"] >= 20,
                        ratio,
                        np.nan,
                    ),
                    "plot_errors": np.where(
                        result[split_kind]["noPrehydro"]["raw_by_bin"] >= 20,
                        ratio_error,
                        np.nan,
                    ),
                }
                for variant in VARIANTS:
                    item = result[split_kind][variant]
                    for bin_index in range(len(tau_edges) - 1):
                        spectra_rows.append(
                            {
                                "variant": variant,
                                "radius": f"{radius:.1f}",
                                "pt_low_exclusive_GeV": f"{pt_low:g}",
                                "pt_high_inclusive_GeV": "inf" if pt_high is None else f"{pt_high:g}",
                                "splitting_selection": split_kind,
                                "bin_index": bin_index,
                                "log10_tau_low": f"{tau_edges[bin_index]:.12e}",
                                "log10_tau_high": f"{tau_edges[bin_index + 1]:.12e}",
                                "log10_tau_center": f"{0.5 * (tau_edges[bin_index] + tau_edges[bin_index + 1]):.12e}",
                                "weighted_yield_per_event": f"{item['values'][bin_index]:.12e}",
                                "weighted_yield_stat_error": f"{item['errors'][bin_index]:.12e}",
                                "selected_jet_normalized_density": f"{item['normalized_values'][bin_index]:.12e}",
                                "selected_jet_normalized_density_stat_error": f"{item['normalized_errors'][bin_index]:.12e}",
                                "differential_cross_section_mb": f"{item['cross_section'][bin_index]:.12e}",
                                "cross_section_stat_error_mb": f"{item['cross_section_errors'][bin_index]:.12e}",
                                "pre_over_no_ratio": (
                                    "1.000000000000e+00"
                                    if variant == "noPrehydro"
                                    else f"{ratio[bin_index]:.12e}"
                                ),
                                "ratio_stat_error": (
                                    "0.000000000000e+00"
                                    if variant == "noPrehydro"
                                    else f"{ratio_error[bin_index]:.12e}"
                                ),
                                "absolute_pre_over_no_ratio": (
                                    "1.000000000000e+00"
                                    if variant == "noPrehydro"
                                    else f"{absolute_ratio[bin_index]:.12e}"
                                ),
                                "absolute_ratio_stat_error": (
                                    "0.000000000000e+00"
                                    if variant == "noPrehydro"
                                    else f"{absolute_ratio_error[bin_index]:.12e}"
                                ),
                                "raw_entries_in_bin": int(item["raw_by_bin"][bin_index]),
                                "raw_entries_total": item["raw"],
                                "underflow_entries": item["underflow"],
                                "overflow_entries": item["overflow"],
                            }
                        )

            spectrum_plot(
                result,
                radius,
                pt_low,
                pt_high,
                tau_edges,
                out_dir / f"{args.prefix}_{tag}_formation_time",
                args.sample_label,
            )

            correlation_data: dict[str, dict[str, Any]] = {}
            for name, (x_edges, field_mode, _) in correlation_axes(radius).items():
                correlation_data[name] = {"x_edges": x_edges}
                weighted_totals: dict[str, np.ndarray] = {}
                raw_totals: dict[str, np.ndarray] = {}
                selected_jet_weight_totals: dict[str, float] = {}
                for variant in VARIANTS:
                    arrays = variants[variant]
                    selection = selections[variant]
                    field = name
                    x, x_event_indices = flatten_splits(arrays, selection, field)
                    tau, tau_event_indices = flatten_splits(arrays, selection, "tauF")
                    if not np.array_equal(x_event_indices, tau_event_indices):
                        raise ValueError(f"{name}/tau split ordering differs")
                    x = transform_correlation(field_mode, x)
                    y = np.log10(tau)
                    finite = np.isfinite(x) & np.isfinite(y)
                    raw, _, _ = np.histogram2d(
                        x[finite], y[finite], bins=(x_edges, tau_edges)
                    )
                    weighted, _, _ = np.histogram2d(
                        x[finite],
                        y[finite],
                        bins=(x_edges, tau_edges),
                        weights=weights[x_event_indices[finite]],
                    )
                    area = np.diff(x_edges)[:, np.newaxis] * np.diff(tau_edges)[np.newaxis, :]
                    selected_event_indices = ak.to_numpy(
                        arrays.eventIndex[selection]
                    ).astype(np.int64)
                    selected_counts = event_count(
                        selected_event_indices, event_count_total
                    )
                    selected_jet_weight = float(np.sum(weights * selected_counts))
                    density = weighted / selected_jet_weight / area
                    event_normalized_density = weighted / weight_sum / area
                    cross_section_density = factor * weighted / area
                    raw_totals[variant] = raw
                    weighted_totals[variant] = weighted
                    selected_jet_weight_totals[variant] = selected_jet_weight
                    correlation_data[name][variant] = {"density": density}
                    for x_index in range(len(x_edges) - 1):
                        for y_index in range(len(tau_edges) - 1):
                            correlation_rows.append(
                                {
                                    "variant": variant,
                                    "radius": f"{radius:.1f}",
                                    "pt_low_exclusive_GeV": f"{pt_low:g}",
                                    "pt_high_inclusive_GeV": "inf" if pt_high is None else f"{pt_high:g}",
                                    "correlation": name,
                                    "x_bin_index": x_index,
                                    "x_low": f"{x_edges[x_index]:.12e}",
                                    "x_high": f"{x_edges[x_index + 1]:.12e}",
                                    "log10_tau_bin_index": y_index,
                                    "log10_tau_low": f"{tau_edges[y_index]:.12e}",
                                    "log10_tau_high": f"{tau_edges[y_index + 1]:.12e}",
                                    "raw_entries": int(raw[x_index, y_index]),
                                    "selected_jet_normalized_density": f"{density[x_index, y_index]:.12e}",
                                    "weighted_density_per_event": f"{event_normalized_density[x_index, y_index]:.12e}",
                                    "weighted_density_mb": f"{cross_section_density[x_index, y_index]:.12e}",
                                }
                            )
                correlation_data[name]["no_raw"] = raw_totals["noPrehydro"]
                correlation_data[name]["ratio"] = np.divide(
                    weighted_totals["withPrehydro"]
                    * selected_jet_weight_totals["noPrehydro"],
                    weighted_totals["noPrehydro"]
                    * selected_jet_weight_totals["withPrehydro"],
                    out=np.full_like(weighted_totals["withPrehydro"], np.nan),
                    where=(weighted_totals["noPrehydro"] > 0.0)
                    & (selected_jet_weight_totals["withPrehydro"] > 0.0),
                )
            correlation_plot(
                correlation_data,
                radius,
                pt_low,
                pt_high,
                tau_edges,
                out_dir / f"{args.prefix}_{tag}_formation_time_correlations",
                args.sample_label,
            )
            metadata["selections"][tag] = selection_metadata

    write_tsv(out_dir / f"{args.prefix}_spectra.tsv", spectra_rows)
    write_tsv(out_dir / f"{args.prefix}_summary.tsv", summary_rows)
    write_tsv(out_dir / f"{args.prefix}_correlations.tsv", correlation_rows)
    metadata_path = out_dir / f"{args.prefix}_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return metadata


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--input-root", required=True, type=Path)
    result.add_argument("--out-dir", required=True, type=Path)
    result.add_argument("--prefix", default="oo5360_jet_formation_time")
    result.add_argument("--sample-label", default="paired OO sample")
    result.add_argument("--log-tau-min", type=float, default=-4.0)
    result.add_argument("--log-tau-max", type=float, default=7.0)
    result.add_argument("--log-tau-bins", type=int, default=44)
    return result


def main() -> int:
    args = parser().parse_args()
    if args.log_tau_max <= args.log_tau_min or args.log_tau_bins <= 0:
        raise ValueError("invalid log-tau binning")
    analyze(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
