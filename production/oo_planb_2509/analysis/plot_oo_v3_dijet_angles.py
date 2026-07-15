#!/usr/bin/env python3
"""Plot same-seed V3 leading/subleading-jet angular correlations."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import awkward as ak
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpecFromSubplotSpec
import numpy as np
import uproot

import plot_oo_jet_variables as base


VARIANTS = ("no_prehydro", "prehydro_alpha037", "prehydro_alpha0335")
RATIO_VARIANTS = ("prehydro_alpha037", "prehydro_alpha0335")
STYLES = {
    "no_prehydro": {
        "label": "No pre-hydro",
        "color": "#111111",
    },
    "prehydro_alpha037": {
        "label": r"Pre-hydro, $\alpha=0.37$",
        "color": "#0072B2",
    },
    "prehydro_alpha0335": {
        "label": r"Pre-hydro, $\alpha=0.335$",
        "color": "#D62728",
    },
}
PAIR_FIELDS = (
    "pairId",
    "sourceIndex",
    "chunkId",
    "hydroIndex",
    "seed",
    "eventNumber",
    "eventWeight",
    "sigmaGen",
    "hardX",
    "hardY",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def exact_numpy_equal(first: np.ndarray, second: np.ndarray) -> bool:
    return first.shape == second.shape and np.array_equal(first, second, equal_nan=True)


def exact_jagged_equal(first: ak.Array, second: ak.Array) -> bool:
    first_counts = ak.to_numpy(ak.num(first, axis=1))
    second_counts = ak.to_numpy(ak.num(second, axis=1))
    if not exact_numpy_equal(first_counts, second_counts):
        return False
    return exact_numpy_equal(
        ak.to_numpy(ak.flatten(first)), ak.to_numpy(ak.flatten(second))
    )


def assert_pair_identity(first: dict[str, np.ndarray], second: dict[str, np.ndarray]) -> None:
    if set(first) != set(second):
        raise ValueError("the two Pairs trees expose different identity fields")
    for field in first:
        if not exact_numpy_equal(first[field], second[field]):
            raise ValueError(f"the two input roots differ in Pairs/{field}")


def weighted_quantile_edges(
    values: np.ndarray,
    weights: np.ndarray,
    *,
    lower: float,
    upper: float,
    bins: int,
) -> np.ndarray:
    """Build common baseline-weighted quantile bins with fixed physical bounds."""
    finite = (
        np.isfinite(values)
        & np.isfinite(weights)
        & (weights > 0.0)
        & (values >= lower)
        & (values <= upper)
    )
    values = values[finite]
    weights = weights[finite]
    if len(values) < bins:
        raise ValueError("too few finite entries for weighted-quantile binning")
    order = np.argsort(values, kind="stable")
    ordered_values = values[order]
    ordered_weights = weights[order]
    cumulative = np.cumsum(ordered_weights)
    cumulative /= cumulative[-1]
    interior = np.interp(
        np.arange(1, bins, dtype=float) / bins,
        cumulative,
        ordered_values,
    )
    edges = np.concatenate(([lower], interior, [upper]))
    if np.any(np.diff(edges) <= 0.0):
        raise ValueError("weighted-quantile bin edges are not strictly increasing")
    return edges


def load_inputs(
    alpha037_root: Path,
    alpha0335_root: Path,
    radius_digits: tuple[int, ...],
) -> tuple[dict[str, np.ndarray], dict[str, ak.Array]]:
    branches = ["pairId", "eventWeight", "sigmaGen"]
    for radius_digit in radius_digits:
        branches.extend(
            (
                f"jet{radius_digit}Pt",
                f"jet{radius_digit}Eta",
                f"jet{radius_digit}Phi",
            )
        )

    with uproot.open(alpha037_root) as first_file:
        first_pairs = first_file["Pairs"].arrays(PAIR_FIELDS, library="np")
        no_prehydro = first_file["noPrehydro/Jets"].arrays(branches, library="ak")
        prehydro_alpha037 = first_file["withPrehydro/Jets"].arrays(
            branches, library="ak"
        )
    with uproot.open(alpha0335_root) as second_file:
        second_pairs = second_file["Pairs"].arrays(PAIR_FIELDS, library="np")
        duplicate_no_prehydro = second_file["noPrehydro/Jets"].arrays(
            branches, library="ak"
        )
        prehydro_alpha0335 = second_file["withPrehydro/Jets"].arrays(
            branches, library="ak"
        )

    assert_pair_identity(first_pairs, second_pairs)
    pair_ids = first_pairs["pairId"]
    if len(pair_ids) == 0 or len(np.unique(pair_ids)) != len(pair_ids):
        raise ValueError("Pairs/pairId is empty or non-unique")
    for field in branches:
        if field in {"pairId", "eventWeight", "sigmaGen"}:
            first_values = ak.to_numpy(no_prehydro[field])
            second_values = ak.to_numpy(duplicate_no_prehydro[field])
            same = exact_numpy_equal(first_values, second_values)
        else:
            same = exact_jagged_equal(
                no_prehydro[field], duplicate_no_prehydro[field]
            )
        if not same:
            raise ValueError(f"no-prehydro duplication failed for {field}")

    for variant, arrays in (
        ("no_prehydro", no_prehydro),
        ("prehydro_alpha037", prehydro_alpha037),
        ("prehydro_alpha0335", prehydro_alpha0335),
    ):
        if not exact_numpy_equal(ak.to_numpy(arrays["pairId"]), pair_ids):
            raise ValueError(f"{variant} pairId ordering differs from Pairs")
        if not exact_numpy_equal(
            ak.to_numpy(arrays["eventWeight"]), first_pairs["eventWeight"]
        ):
            raise ValueError(f"{variant} eventWeight differs from Pairs")
        if not exact_numpy_equal(
            ak.to_numpy(arrays["sigmaGen"]), first_pairs["sigmaGen"]
        ):
            raise ValueError(f"{variant} sigmaGen differs from Pairs")

    return first_pairs, {
        "no_prehydro": no_prehydro,
        "prehydro_alpha037": prehydro_alpha037,
        "prehydro_alpha0335": prehydro_alpha0335,
    }


def extract_dijet(
    arrays: ak.Array,
    radius_digit: int,
    lead_pt_min: float,
    sublead_pt_min: float,
) -> dict[str, np.ndarray]:
    pt = arrays[f"jet{radius_digit}Pt"]
    order = ak.argsort(pt, axis=1, ascending=False)
    sorted_pt = ak.pad_none(pt[order], 2, axis=1, clip=True)
    sorted_eta = ak.pad_none(
        arrays[f"jet{radius_digit}Eta"][order], 2, axis=1, clip=True
    )
    sorted_phi = ak.pad_none(
        arrays[f"jet{radius_digit}Phi"][order], 2, axis=1, clip=True
    )
    selected = ak.fill_none(
        (sorted_pt[:, 0] > lead_pt_min)
        & (sorted_pt[:, 1] > sublead_pt_min)
        & np.isfinite(sorted_pt[:, 0])
        & np.isfinite(sorted_pt[:, 1])
        & np.isfinite(sorted_eta[:, 0])
        & np.isfinite(sorted_eta[:, 1])
        & np.isfinite(sorted_phi[:, 0])
        & np.isfinite(sorted_phi[:, 1]),
        False,
    )
    selected_numpy = ak.to_numpy(selected).astype(bool)
    event_indices = np.flatnonzero(selected_numpy)
    eta_lead = ak.to_numpy(sorted_eta[selected, 0])
    eta_sublead = ak.to_numpy(sorted_eta[selected, 1])
    phi_lead = ak.to_numpy(sorted_phi[selected, 0])
    phi_sublead = ak.to_numpy(sorted_phi[selected, 1])
    delta_phi = np.abs(
        np.arctan2(
            np.sin(phi_lead - phi_sublead),
            np.cos(phi_lead - phi_sublead),
        )
    )
    return {
        "selected": selected_numpy.astype(np.float64),
        "event_indices": event_indices,
        "abs_delta_eta": np.abs(eta_lead - eta_sublead),
        "delta_phi": delta_phi,
    }


def normalized_histogram(
    values: np.ndarray,
    event_indices: np.ndarray,
    selected: np.ndarray,
    weights: np.ndarray,
    edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    matrix, underflow, overflow = base.histogram_matrix(
        values, event_indices, weights, edges
    )
    selected_by_event = weights * selected
    selected_total = float(np.sum(selected_by_event))
    if selected_total <= 0.0:
        raise ValueError("dijet selection has zero weighted entries")
    widths = np.diff(edges)
    totals = np.sum(matrix, axis=0)
    density = totals / selected_total / widths
    leave_selected = selected_total - selected_by_event
    leave_density = np.divide(
        totals[np.newaxis, :] - matrix,
        leave_selected[:, np.newaxis],
        out=np.full_like(matrix, np.nan),
        where=leave_selected[:, np.newaxis] > 0.0,
    ) / widths[np.newaxis, :]
    return density, base.jackknife_error(leave_density), matrix, underflow, overflow


def plot_panel(
    figure,
    outer_spec,
    *,
    radius: float,
    variable: str,
    label: str,
    edges: np.ndarray,
    results: dict[str, dict[str, np.ndarray]],
    show_legend: bool,
) -> None:
    inner = GridSpecFromSubplotSpec(
        2,
        1,
        subplot_spec=outer_spec,
        height_ratios=(3.0, 1.0),
        hspace=0.05,
    )
    upper = figure.add_subplot(inner[0])
    lower = figure.add_subplot(inner[1], sharex=upper)
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    for variant in VARIANTS:
        style = STYLES[variant]
        values = results[variant][f"{variable}_density"]
        errors = results[variant][f"{variable}_error"]
        finite = np.isfinite(values) & np.isfinite(errors)
        if variant == "no_prehydro":
            upper.stairs(
                values,
                edges,
                color=style["color"],
                linewidth=1.6,
                label=style["label"],
            )
        else:
            shift = -0.07 if variant == "prehydro_alpha037" else 0.07
            upper.errorbar(
                centers[finite] + shift * widths[finite],
                values[finite],
                yerr=errors[finite],
                color=style["color"],
                marker="o",
                markersize=3.1,
                linestyle="none",
                linewidth=0.8,
                capsize=1.4,
                label=style["label"],
            )
    upper.set_ylabel(r"$(1/\sigma_{\rm dijet})\,d\sigma/dx$", fontsize=8.2)
    upper.tick_params(labelbottom=False, labelsize=7.7)
    upper.grid(alpha=0.22)
    upper.text(
        0.03,
        0.92,
        rf"anti-$k_T$ $R={radius:.1f}$",
        transform=upper.transAxes,
        ha="left",
        va="top",
        fontsize=8.5,
    )
    if show_legend:
        upper.legend(frameon=False, fontsize=7.5, loc="best")

    lower.axhline(1.0, color="0.45", linewidth=0.9)
    all_ratios = []
    all_errors = []
    for variant in RATIO_VARIANTS:
        style = STYLES[variant]
        ratio = results[variant][f"{variable}_ratio"]
        error = results[variant][f"{variable}_ratio_error"]
        all_ratios.append(ratio)
        all_errors.append(error)
        finite = np.isfinite(ratio) & np.isfinite(error)
        shift = -0.07 if variant == "prehydro_alpha037" else 0.07
        lower.errorbar(
            centers[finite] + shift * widths[finite],
            ratio[finite],
            yerr=error[finite],
            color=style["color"],
            marker="o",
            markersize=2.7,
            linestyle="none",
            linewidth=0.75,
            capsize=1.2,
        )
    lower.set_ylim(
        *base.ratio_limits(np.concatenate(all_ratios), np.concatenate(all_errors))
    )
    lower.set_ylabel("Pre/no", fontsize=7.5)
    lower.set_xlabel(label, fontsize=8.6)
    lower.tick_params(labelsize=7.5)
    lower.grid(alpha=0.22)


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def run(args: argparse.Namespace) -> dict[str, object]:
    if args.lead_pt_min <= args.sublead_pt_min:
        raise ValueError("--lead-pt-min must exceed --sublead-pt-min")
    if args.sublead_pt_min <= 0.0:
        raise ValueError("--sublead-pt-min must be positive")
    if args.bins < 2:
        raise ValueError("--bins must be at least two")
    alpha037_root = args.alpha037_root.expanduser().resolve()
    alpha0335_root = args.alpha0335_root.expanduser().resolve()
    out_dir = args.out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    radius_digits = tuple(args.radius_digits)
    if not radius_digits or any(radius <= 0 for radius in radius_digits):
        raise ValueError("--radius-digits must contain positive integers")

    pairs, arrays = load_inputs(alpha037_root, alpha0335_root, radius_digits)
    weights = pairs["eventWeight"].astype(np.float64)
    sigma_gen = pairs["sigmaGen"].astype(np.float64)
    if not np.all(np.isfinite(weights)) or not np.all(weights > 0.0):
        raise ValueError("event weights must be finite and positive")
    if not np.all(np.isfinite(sigma_gen)) or not np.all(sigma_gen > 0.0):
        raise ValueError("sigmaGen must be finite and positive")
    _, weight_sum, sigma_merged = base.pythia_parallel_factor(weights, sigma_gen)

    all_results: dict[int, dict[str, dict[str, np.ndarray]]] = {}
    histogram_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    metadata_radii: dict[str, object] = {}
    variable_definitions = {
        "abs_delta_eta": {
            "label": r"$|\Delta\eta_{12}|$",
            "lower": 0.0,
            "upper": 10.0,
        },
        "delta_phi": {
            "label": r"$|\Delta\phi_{12}|$ [rad]",
            "lower": 0.0,
            "upper": math.pi,
        },
    }

    for radius_digit in radius_digits:
        extracted = {
            variant: extract_dijet(
                arrays[variant],
                radius_digit,
                args.lead_pt_min,
                args.sublead_pt_min,
            )
            for variant in VARIANTS
        }
        edges = {
            variable: weighted_quantile_edges(
                extracted["no_prehydro"][variable],
                weights[extracted["no_prehydro"]["event_indices"]],
                lower=definition["lower"],
                upper=definition["upper"],
                bins=args.bins,
            )
            for variable, definition in variable_definitions.items()
        }
        radius_results: dict[str, dict[str, np.ndarray]] = {
            variant: {} for variant in VARIANTS
        }
        matrices: dict[str, dict[str, np.ndarray]] = {
            variant: {} for variant in VARIANTS
        }
        flow: dict[str, dict[str, tuple[int, int]]] = {
            variant: {} for variant in VARIANTS
        }
        for variant in VARIANTS:
            for variable in variable_definitions:
                density, error, matrix, underflow, overflow = normalized_histogram(
                    extracted[variant][variable],
                    extracted[variant]["event_indices"],
                    extracted[variant]["selected"],
                    weights,
                    edges[variable],
                )
                radius_results[variant][f"{variable}_density"] = density
                radius_results[variant][f"{variable}_error"] = error
                matrices[variant][variable] = matrix
                flow[variant][variable] = (underflow, overflow)

        for variable in variable_definitions:
            for variant in RATIO_VARIANTS:
                ratio, ratio_error = base.paired_normalized_ratio(
                    matrices[variant][variable],
                    matrices["no_prehydro"][variable],
                    extracted[variant]["selected"],
                    extracted["no_prehydro"]["selected"],
                    weights,
                )
                radius_results[variant][f"{variable}_ratio"] = ratio
                radius_results[variant][f"{variable}_ratio_error"] = ratio_error
            radius_results["no_prehydro"][f"{variable}_ratio"] = np.ones(
                args.bins
            )
            radius_results["no_prehydro"][f"{variable}_ratio_error"] = np.zeros(
                args.bins
            )

            for variant in VARIANTS:
                for bin_index in range(args.bins):
                    histogram_rows.append(
                        {
                            "radius": f"{radius_digit / 10.0:.1f}",
                            "variable": variable,
                            "variant": variant,
                            "bin_index": bin_index,
                            "bin_low": f"{edges[variable][bin_index]:.12e}",
                            "bin_high": f"{edges[variable][bin_index + 1]:.12e}",
                            "normalized_density": f"{radius_results[variant][f'{variable}_density'][bin_index]:.12e}",
                            "normalized_density_stat_error": f"{radius_results[variant][f'{variable}_error'][bin_index]:.12e}",
                            "ratio_to_no_prehydro": f"{radius_results[variant][f'{variable}_ratio'][bin_index]:.12e}",
                            "ratio_stat_error": f"{radius_results[variant][f'{variable}_ratio_error'][bin_index]:.12e}",
                        }
                    )

        radius_summary: dict[str, object] = {
            "weightedQuantileEdgesFromNoPrehydro": {
                variable: edge_values.tolist()
                for variable, edge_values in edges.items()
            },
            "variants": {},
        }
        for variant in VARIANTS:
            selected = extracted[variant]["selected"]
            selected_weight = float(np.sum(weights * selected))
            if variant == "no_prehydro":
                yield_ratio, yield_ratio_error = 1.0, 0.0
            else:
                yield_ratio, yield_ratio_error = base.paired_integrated_ratio(
                    selected,
                    extracted["no_prehydro"]["selected"],
                    weights,
                )
            summary_rows.append(
                {
                    "radius": f"{radius_digit / 10.0:.1f}",
                    "variant": variant,
                    "raw_selected_dijets": int(np.sum(selected)),
                    "weighted_selected_dijets": f"{selected_weight:.12e}",
                    "weighted_selected_fraction": f"{selected_weight / weight_sum:.12e}",
                    "selected_yield_ratio_to_no_prehydro": f"{yield_ratio:.12e}",
                    "selected_yield_ratio_stat_error": f"{yield_ratio_error:.12e}",
                }
            )
            radius_summary["variants"][variant] = {
                "rawSelectedDijets": int(np.sum(selected)),
                "weightedSelectedDijets": selected_weight,
                "weightedSelectedFraction": selected_weight / weight_sum,
                "selectedYieldRatioToNoPrehydro": yield_ratio,
                "selectedYieldRatioError": yield_ratio_error,
                "histogramFlow": {
                    variable: {
                        "underflow": flow[variant][variable][0],
                        "overflow": flow[variant][variable][1],
                    }
                    for variable in variable_definitions
                },
            }
        metadata_radii[f"{radius_digit / 10.0:.1f}"] = radius_summary
        all_results[radius_digit] = radius_results

    figure = plt.figure(figsize=(11.6, 8.2))
    grid = figure.add_gridspec(
        len(radius_digits), 2, wspace=0.29, hspace=0.30
    )
    for row_index, radius_digit in enumerate(radius_digits):
        radius = radius_digit / 10.0
        for column_index, (variable, definition) in enumerate(
            variable_definitions.items()
        ):
            plot_panel(
                figure,
                grid[row_index, column_index],
                radius=radius,
                variable=variable,
                label=definition["label"],
                edges=np.asarray(
                    metadata_radii[f"{radius:.1f}"][
                        "weightedQuantileEdgesFromNoPrehydro"
                    ][variable]
                ),
                results=all_results[radius_digit],
                show_legend=row_index == 0 and column_index == 0,
            )
    figure.suptitle(
        r"O+O 5.36 TeV, 0--5%, V3 same-seed leading/subleading-jet angles",
        fontsize=13.0,
        y=0.992,
    )
    figure.text(
        0.5,
        0.012,
        rf"4MomSub jets; $p_{{T,1}}>{args.lead_pt_min:g}$ GeV, "
        rf"$p_{{T,2}}>{args.sublead_pt_min:g}$ GeV; no added $\eta$ cut. "
        r"Per-variant area normalization; baseline-weighted decile bins; paired delete-one-event jackknife.",
        ha="center",
        fontsize=8.2,
    )
    figure.subplots_adjust(top=0.94, bottom=0.075)

    plot_pdf = out_dir / f"{args.prefix}.pdf"
    plot_png = plot_pdf.with_suffix(".png")
    figure.savefig(plot_pdf)
    figure.savefig(plot_png, dpi=180)
    plt.close(figure)

    radius_plots: dict[str, dict[str, str]] = {}
    for radius_digit in radius_digits:
        radius = radius_digit / 10.0
        radius_figure = plt.figure(figsize=(11.8, 4.7))
        radius_grid = radius_figure.add_gridspec(1, 2, wspace=0.29)
        for column_index, (variable, definition) in enumerate(
            variable_definitions.items()
        ):
            plot_panel(
                radius_figure,
                radius_grid[0, column_index],
                radius=radius,
                variable=variable,
                label=definition["label"],
                edges=np.asarray(
                    metadata_radii[f"{radius:.1f}"][
                        "weightedQuantileEdgesFromNoPrehydro"
                    ][variable]
                ),
                results=all_results[radius_digit],
                show_legend=column_index == 0,
            )
        radius_figure.suptitle(
            rf"O+O 5.36 TeV, 0--5%, V3 same-seed dijet angles; anti-$k_T$ $R={radius:.1f}$",
            fontsize=12.7,
            y=0.992,
        )
        radius_figure.text(
            0.5,
            0.015,
            rf"4MomSub jets; $p_{{T,1}}>{args.lead_pt_min:g}$ GeV, "
            rf"$p_{{T,2}}>{args.sublead_pt_min:g}$ GeV; no added $\eta$ cut. "
            r"Per-variant area normalization; baseline-weighted decile bins; paired delete-one-event jackknife.",
            ha="center",
            fontsize=8.1,
        )
        radius_figure.subplots_adjust(top=0.90, bottom=0.14)
        radius_pdf = out_dir / f"{args.prefix}_R0{radius_digit}.pdf"
        radius_png = radius_pdf.with_suffix(".png")
        radius_figure.savefig(radius_pdf)
        radius_figure.savefig(radius_png, dpi=180)
        plt.close(radius_figure)
        radius_plots[f"{radius:.1f}"] = {
            "pdf": str(radius_pdf),
            "png": str(radius_png),
        }

    histogram_path = out_dir / f"{args.prefix}_histograms.tsv"
    summary_path = out_dir / f"{args.prefix}_summary.tsv"
    write_tsv(histogram_path, histogram_rows)
    write_tsv(summary_path, summary_rows)
    metadata = {
        "status": "PASS",
        "publicationStatus": "PROVISIONAL_FIXED_BOUNDARY_SAME_SEED_DIAGNOSTIC",
        "pairCount": len(pairs["pairId"]),
        "uniqueSeedCount": len(np.unique(pairs["seed"])),
        "variants": list(VARIANTS),
        "radii": metadata_radii,
        "selection": {
            "jetDefinition": "anti-kT, 4MomSub-corrected four-vector",
            "ordering": "descending corrected jet pT independently in each variant",
            "leadPtMinGeVExclusive": args.lead_pt_min,
            "subleadPtMinGeVExclusive": args.sublead_pt_min,
            "additionalJetEtaCut": None,
            "absDeltaEta": "abs(eta_lead - eta_sublead) using corrected jet eta",
            "absDeltaPhi": "abs(atan2(sin(phi_lead-phi_sublead), cos(phi_lead-phi_sublead))) in [0,pi]",
        },
        "normalization": (
            "eventWeight is applied exactly once; each variant is normalized by "
            "its own weighted selected-dijet yield, so every displayed angular "
            "density has unit area"
        ),
        "binning": (
            "ten weighted-quantile bins derived separately for each radius and "
            "observable from the no-prehydro selected dijets; the same fixed edges "
            "are then used for all three variants"
        ),
        "shapeRatio": (
            "area-normalized prehydro angular density divided by the area-normalized "
            "no-prehydro density; each normalization is recomputed in every paired replica"
        ),
        "uncertainty": "paired delete-one-event jackknife over the aligned hard-event rows",
        "weightSum": weight_sum,
        "sigmaMergedMb": sigma_merged,
        "inputs": {
            "alpha037Root": str(alpha037_root),
            "alpha037RootSha256": sha256(alpha037_root),
            "alpha0335Root": str(alpha0335_root),
            "alpha0335RootSha256": sha256(alpha0335_root),
        },
        "outputs": {
            "plotPdf": str(plot_pdf),
            "plotPng": str(plot_png),
            "radiusPlots": radius_plots,
            "histograms": str(histogram_path),
            "summary": str(summary_path),
        },
    }
    metadata_path = out_dir / f"{args.prefix}_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    print(json.dumps(metadata, indent=2, sort_keys=True))
    return metadata


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--alpha037-root", type=Path, required=True)
    parser.add_argument("--alpha0335-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--lead-pt-min", type=float, default=30.0)
    parser.add_argument("--sublead-pt-min", type=float, default=20.0)
    parser.add_argument("--radius-digits", type=int, nargs="+", default=(4, 8))
    parser.add_argument("--bins", type=int, default=10)
    parser.add_argument("--prefix", default="oo5360_v3_dijet_angles")
    return parser


def main() -> int:
    run(make_parser().parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
