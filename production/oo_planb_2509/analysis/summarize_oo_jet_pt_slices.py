#!/usr/bin/env python3
"""Summarize and validate bounded-pT OO jet comparisons."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


VARIANTS = ("noPrehydro", "withPrehydro")
VARIANT_LABELS = {
    "noPrehydro": "No pre-hydro",
    "withPrehydro": "Plan B pre-hydro",
}
VARIANT_COLORS = {
    "noPrehydro": "#0072B2",
    "withPrehydro": "#D55E00",
}


def load_metadata(path: Path) -> dict[str, object]:
    return json.loads(path.read_text())


def range_label(metadata: dict[str, object]) -> str:
    low = float(metadata["ptMinGeV"])
    high = metadata["ptMaxGeV"]
    if high is None:
        return f">{low:g}"
    return f"{low:g}-{float(high):g}"


def validate_slices(
    slices: list[dict[str, object]], inclusive: dict[str, object]
) -> dict[str, object]:
    if not slices:
        raise ValueError("at least one pT slice is required")
    slices.sort(key=lambda item: float(item["ptMinGeV"]))
    checks: list[str] = []

    reference_fields = (
        "inputRoot",
        "inputBytes",
        "pairCount",
        "uniqueSeedCount",
        "weightSum",
        "weightedSigmaSum",
        "sigmaMergedMb",
        "crossSectionFactor",
        "effectiveEventCountFromWeights",
        "normalization",
        "uncertainty",
    )
    for field in reference_fields:
        expected = inclusive[field]
        if any(item[field] != expected for item in slices):
            raise ValueError(f"slice metadata disagree on {field}")
    checks.append("all slices use the inclusive sample and normalization")

    if float(slices[0]["ptMinGeV"]) != float(inclusive["ptMinGeV"]):
        raise ValueError("first slice does not start at the inclusive threshold")
    for lower, upper in zip(slices[:-1], slices[1:]):
        if lower["ptMaxGeV"] is None or float(lower["ptMaxGeV"]) != float(upper["ptMinGeV"]):
            raise ValueError("pT slices have a gap or overlap")
    if slices[-1]["ptMaxGeV"] is not None:
        raise ValueError("last pT slice is not open ended")
    checks.append("slice boundaries form one exclusive-lower/inclusive-upper partition")

    closure: dict[str, dict[str, float]] = {}
    for radius in ("0.4", "0.8"):
        closure[radius] = {}
        for variant in VARIANTS:
            sliced = sum(
                float(item["radii"][radius]["variants"][variant]["integratedJetCrossSectionMb"])
                for item in slices
            )
            expected = float(
                inclusive["radii"][radius]["variants"][variant]["integratedJetCrossSectionMb"]
            )
            difference = sliced - expected
            if not math.isclose(sliced, expected, rel_tol=1.0e-12, abs_tol=1.0e-14):
                raise ValueError(f"R={radius} {variant} slices do not close to inclusive")
            closure[radius][variant] = difference
    checks.append("all four sliced cross-section sums reproduce the inclusive result")

    return {
        "status": "PASS",
        "checks": checks,
        "crossSectionClosureDifferenceMb": closure,
    }


def write_summary(
    path: Path, slices: list[dict[str, object]], validation: dict[str, object]
) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "pt_range_gev",
                "pt_min_gev_exclusive",
                "pt_max_gev_inclusive",
                "radius",
                "variant",
                "integrated_jet_cross_section_mb",
                "stat_error_mb",
                "pre_over_no_ratio",
                "ratio_stat_error",
                "raw_selected_jets",
                "weighted_selected_jets",
            ]
        )
        for item in slices:
            for radius in ("0.4", "0.8"):
                ratio = item["radii"][radius]["integratedPreOverNoRatio"]
                ratio_error = item["radii"][radius]["integratedPreOverNoRatioError"]
                for variant in VARIANTS:
                    result = item["radii"][radius]["variants"][variant]
                    writer.writerow(
                        [
                            range_label(item),
                            item["ptMinGeV"],
                            "inf" if item["ptMaxGeV"] is None else item["ptMaxGeV"],
                            radius,
                            variant,
                            f'{result["integratedJetCrossSectionMb"]:.12e}',
                            f'{result["integratedJetCrossSectionErrorMb"]:.12e}',
                            f"{ratio:.12e}",
                            f"{ratio_error:.12e}",
                            result["rawSelectedJets"],
                            f'{result["weightedSelectedJets"]:.12e}',
                        ]
                    )
    validation["summaryTsv"] = path.name


def plot_summary(out_path: Path, slices: list[dict[str, object]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [range_label(item) for item in slices]
    x = np.arange(len(slices), dtype=float)
    figure, axes = plt.subplots(
        2,
        2,
        figsize=(11.5, 7.6),
        sharex="col",
        gridspec_kw={"height_ratios": [2.25, 1.0], "hspace": 0.08, "wspace": 0.24},
    )
    for column, radius in enumerate(("0.4", "0.8")):
        upper = axes[0, column]
        lower = axes[1, column]
        for offset, variant in zip((-0.08, 0.08), VARIANTS):
            values = [
                item["radii"][radius]["variants"][variant]["integratedJetCrossSectionMb"]
                for item in slices
            ]
            errors = [
                item["radii"][radius]["variants"][variant]["integratedJetCrossSectionErrorMb"]
                for item in slices
            ]
            upper.errorbar(
                x + offset,
                values,
                yerr=errors,
                color=VARIANT_COLORS[variant],
                marker="o" if variant == "noPrehydro" else "s",
                markersize=5.0,
                linewidth=1.2,
                capsize=2.5,
                label=VARIANT_LABELS[variant],
            )
        ratios = [item["radii"][radius]["integratedPreOverNoRatio"] for item in slices]
        ratio_errors = [
            item["radii"][radius]["integratedPreOverNoRatioError"] for item in slices
        ]
        lower.axhline(1.0, color="0.45", linewidth=1.0)
        lower.errorbar(
            x,
            ratios,
            yerr=ratio_errors,
            color="#333333",
            marker="o",
            markersize=4.8,
            linestyle="none",
            capsize=2.5,
        )
        upper.set_yscale("log")
        upper.set_title(rf"anti-$k_T$ $R={radius}$", fontsize=12)
        upper.grid(alpha=0.22)
        lower.grid(alpha=0.22)
        lower.set_ylim(0.92, 1.02)
        lower.set_xticks(x, labels)
        lower.set_xlabel(r"corrected jet $p_T$ interval [GeV]")
        lower.set_ylabel("Plan B / no pre-hydro")
        if column == 0:
            upper.set_ylabel(r"integrated jet cross section [mb]")
            upper.legend(frameon=False, fontsize=9)
    figure.suptitle(
        r"O16+O16 5.36 TeV, 0--5% diagnostic: weighted jet yields by $p_T$ interval",
        fontsize=14,
        y=0.985,
    )
    figure.text(
        0.5,
        0.012,
        "PythiaParallel merged normalization; paired delete-one-run jackknife; no additional jet-eta cut.",
        ha="center",
        fontsize=8.5,
    )
    figure.subplots_adjust(top=0.91, bottom=0.12)
    figure.savefig(out_path.with_suffix(".pdf"))
    figure.savefig(out_path.with_suffix(".png"), dpi=180)
    plt.close(figure)


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slice-metadata", nargs="+", type=Path, required=True)
    parser.add_argument("--inclusive-metadata", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="oo5360_jet_pt_slice_summary")
    return parser


def main() -> int:
    args = make_parser().parse_args()
    slices = [load_metadata(path) for path in args.slice_metadata]
    inclusive = load_metadata(args.inclusive_metadata)
    validation = validate_slices(slices, inclusive)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem = args.out_dir / args.prefix
    write_summary(stem.with_suffix(".tsv"), slices, validation)
    plot_summary(stem, slices)
    validation["sliceMetadata"] = [path.name for path in args.slice_metadata]
    validation["inclusiveMetadata"] = args.inclusive_metadata.name
    validation_path = stem.with_name(f"{stem.name}_validation.json")
    validation_path.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
