#!/usr/bin/env python3
"""Validate and summarize strict three-way V3 OO jet-pT slices."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


VARIANTS = ("no_prehydro", "prehydro_alpha037", "prehydro_alpha0335")
RATIO_VARIANTS = ("prehydro_alpha037", "prehydro_alpha0335")
RADII = ("0.1", "0.2", "0.4", "0.8")
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


def load(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text())
    if value.get("status") != "PASS":
        raise ValueError(f"input metadata is not PASS: {path}")
    return value


def range_label(metadata: dict[str, object]) -> str:
    low = float(metadata["ptMinGeV"])
    high = metadata["ptMaxGeV"]
    if high is None:
        return f">{low:g}"
    return f"{low:g}-{float(high):g}"


def value(
    metadata: dict[str, object], radius: str, variant: str, field: str
) -> float:
    return float(metadata["radii"][radius]["variants"][variant][field])


def validate(
    slices: list[dict[str, object]], inclusive: dict[str, object]
) -> dict[str, object]:
    if len(slices) != 4:
        raise ValueError(f"expected four disjoint pT slices, found {len(slices)}")
    slices.sort(key=lambda item: float(item["ptMinGeV"]))
    common_fields = (
        "publicationStatus",
        "expectedEvents",
        "normalization",
        "uncertainty",
    )
    for field in common_fields:
        if any(item[field] != inclusive[field] for item in slices):
            raise ValueError(f"slice metadata disagree on {field}")
    input_hashes = {
        key: value
        for key, value in inclusive["inputs"].items()
        if key.endswith("Sha256")
    }
    for item in slices:
        candidate = {
            key: value
            for key, value in item["inputs"].items()
            if key.endswith("Sha256")
        }
        # Slice histograms differ, so input hashes are recorded but cannot match.
        if set(candidate) != set(input_hashes):
            raise ValueError("slice metadata have inconsistent input roles")

    if not math.isclose(float(inclusive["ptMinGeV"]), 20.0):
        raise ValueError("inclusive spectrum must start at 20 GeV")
    if inclusive["ptMaxGeV"] is not None:
        raise ValueError("inclusive spectrum must have an open upper boundary")
    if not math.isclose(float(slices[0]["ptMinGeV"]), 20.0):
        raise ValueError("first slice must start at 20 GeV")
    for lower, upper in zip(slices[:-1], slices[1:]):
        if lower["ptMaxGeV"] is None or not math.isclose(
            float(lower["ptMaxGeV"]), float(upper["ptMinGeV"])
        ):
            raise ValueError("pT slices have a gap or overlap")
    if slices[-1]["ptMaxGeV"] is not None:
        raise ValueError("last pT slice must be open ended")

    closure: dict[str, dict[str, float]] = {}
    for radius in RADII:
        closure[radius] = {}
        for variant in VARIANTS:
            sliced = sum(
                value(item, radius, variant, "integratedJetCrossSectionMb")
                for item in slices
            )
            expected = value(
                inclusive, radius, variant, "integratedJetCrossSectionMb"
            )
            difference = sliced - expected
            if not math.isclose(
                sliced, expected, rel_tol=1.0e-11, abs_tol=1.0e-14
            ):
                raise ValueError(
                    f"R={radius} {variant} slices do not close: {sliced} != {expected}"
                )
            closure[radius][variant] = difference
    return {
        "status": "PASS",
        "checks": [
            "four slices form the exclusive-lower/inclusive-upper pT>20 partition",
            "all slices use the same event count and absolute cross-section normalization",
            "all radii and all three variants close to their inclusive spectra",
        ],
        "crossSectionClosureDifferenceMb": closure,
    }


def write_tsv(path: Path, slices: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "pt_range_gev",
                "pt_min_gev_exclusive",
                "pt_max_gev_inclusive",
                "radius",
                "variant",
                "integrated_jet_cross_section_mb",
                "stat_error_mb",
                "ratio_to_no_prehydro",
                "ratio_stat_error",
            ]
        )
        for item in slices:
            for radius in RADII:
                for variant in VARIANTS:
                    writer.writerow(
                        [
                            range_label(item),
                            item["ptMinGeV"],
                            "inf" if item["ptMaxGeV"] is None else item["ptMaxGeV"],
                            radius,
                            variant,
                            f"{value(item, radius, variant, 'integratedJetCrossSectionMb'):.12e}",
                            f"{value(item, radius, variant, 'integratedJetCrossSectionErrorMb'):.12e}",
                            f"{value(item, radius, variant, 'ratioToNoPrehydro'):.12e}",
                            f"{value(item, radius, variant, 'ratioToNoPrehydroError'):.12e}",
                        ]
                    )


def plot(path: Path, slices: list[dict[str, object]]) -> None:
    labels = [range_label(item) for item in slices]
    x = np.arange(len(slices), dtype=float)
    figure, axes = plt.subplots(
        2,
        len(RADII),
        figsize=(20.0, 7.6),
        sharex="col",
        gridspec_kw={"height_ratios": [2.25, 1.0], "hspace": 0.08, "wspace": 0.24},
    )
    for column, radius in enumerate(RADII):
        upper = axes[0, column]
        lower = axes[1, column]
        for offset, variant in zip((-0.12, 0.0, 0.12), VARIANTS):
            style = STYLES[variant]
            values = [
                value(item, radius, variant, "integratedJetCrossSectionMb")
                for item in slices
            ]
            if variant == "no_prehydro":
                upper.plot(
                    x,
                    values,
                    color=style["color"],
                    drawstyle="steps-mid",
                    linewidth=1.7,
                    label=style["label"],
                )
            else:
                upper.errorbar(
                    x + offset,
                    values,
                    yerr=[
                        value(
                            item,
                            radius,
                            variant,
                            "integratedJetCrossSectionErrorMb",
                        )
                        for item in slices
                    ],
                    color=style["color"],
                    marker=style["marker"],
                    linestyle="none",
                    markersize=4.8,
                    linewidth=1.2,
                    capsize=2.3,
                    label=style["label"],
                )
        lower.axhline(1.0, color="0.45", linewidth=0.9)
        for offset, variant in zip((-0.04, 0.04), RATIO_VARIANTS):
            style = STYLES[variant]
            lower.errorbar(
                x + offset,
                [value(item, radius, variant, "ratioToNoPrehydro") for item in slices],
                yerr=[
                    value(item, radius, variant, "ratioToNoPrehydroError")
                    for item in slices
                ],
                color=style["color"],
                marker=style["marker"],
                markersize=4.2,
                linestyle="none",
                capsize=2.1,
            )
        upper.set_yscale("log")
        upper.set_title(rf"anti-$k_T$ $R={radius}$", fontsize=12)
        upper.grid(alpha=0.22)
        lower.grid(alpha=0.22)
        lower.set_ylim(0.80, 1.20)
        lower.set_xticks(x, labels)
        lower.set_xlabel(r"corrected jet $p_T$ interval [GeV]")
        lower.set_ylabel("pre/no")
        if column == 0:
            upper.set_ylabel("integrated jet cross section [mb]")
            upper.legend(frameon=False, fontsize=8.3)
    figure.suptitle(
        r"O+O 5.36 TeV, 0--5%, V3 same-seed matched: weighted jet yields by $p_T$ interval",
        fontsize=14,
        y=0.985,
    )
    figure.text(
        0.5,
        0.012,
        f"{int(slices[0]['expectedEvents']):,} matched tasks; PythiaParallel normalization; paired delete-one-event jackknife.",
        ha="center",
        fontsize=8.5,
    )
    figure.subplots_adjust(top=0.91, bottom=0.12)
    figure.savefig(path)
    figure.savefig(path.with_suffix(".png"), dpi=180)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--slice-metadata", type=Path, nargs=4, required=True)
    parser.add_argument("--inclusive-metadata", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    slices = [load(path) for path in args.slice_metadata]
    inclusive = load(args.inclusive_metadata)
    validation = validate(slices, inclusive)
    table = args.out_dir / f"{args.prefix}.tsv"
    figure = args.out_dir / f"{args.prefix}.pdf"
    write_tsv(table, slices)
    plot(figure, slices)
    validation.update(
        {
            "publicationStatus": inclusive["publicationStatus"],
            "expectedEvents": inclusive["expectedEvents"],
            "inclusiveMetadata": str(args.inclusive_metadata.resolve()),
            "sliceMetadata": [str(path.resolve()) for path in args.slice_metadata],
            "summaryTsv": str(table),
            "plot": str(figure),
        }
    )
    metadata = args.out_dir / f"{args.prefix}_metadata.json"
    metadata.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n")
    print(json.dumps(validation, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
