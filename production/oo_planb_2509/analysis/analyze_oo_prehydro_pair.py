#!/usr/bin/env python3
"""Analyze paired OO no-prehydro/prehydro AA outputs against a common pp sample."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
from pathlib import Path, PurePosixPath
import re
import subprocess
import tarfile
from typing import BinaryIO


DEFAULT_BINS = [4, 5, 6, 8, 10, 15, 20, 30, 40, 60, 80, 110]
CHARGED_ABS = {211, 321, 2212}
AA_VARIANTS = ["no_prehydro", "with_prehydro"]


@dataclass(frozen=True)
class PythiaRun:
    histogram: list[float]
    sigma_gen: float
    weight_sum: float
    event_count: int


@dataclass(frozen=True)
class SpectrumStats:
    values: list[float]
    standard_errors: list[float]
    event_count: int
    run_count: int
    sigma_gen: float
    weight_sum: float


class PythiaAggregate:
    """Merge independent PYTHIA runs using the PythiaParallel convention."""

    def __init__(self, n_bins: int) -> None:
        self.runs: list[PythiaRun] = []
        self.histogram = [0.0] * n_bins
        self.weight_sum = 0.0
        self.weighted_sigma_sum = 0.0
        self.event_count = 0

    def add(self, run: PythiaRun) -> None:
        if len(run.histogram) != len(self.histogram):
            raise ValueError("PYTHIA run has the wrong number of histogram bins")
        self.runs.append(run)
        for i, value in enumerate(run.histogram):
            self.histogram[i] += value
        self.weight_sum += run.weight_sum
        self.weighted_sigma_sum += run.weight_sum * run.sigma_gen
        self.event_count += run.event_count

    @staticmethod
    def _normalized(
        histogram: list[float], weight_sum: float, weighted_sigma_sum: float
    ) -> tuple[list[float], float]:
        if weight_sum == 0.0:
            return [math.nan] * len(histogram), math.nan
        sigma_gen = weighted_sigma_sum / weight_sum
        scale = sigma_gen / weight_sum
        return [scale * value for value in histogram], sigma_gen

    def stats(self) -> SpectrumStats:
        values, sigma_gen = self._normalized(
            self.histogram, self.weight_sum, self.weighted_sigma_sum
        )
        run_count = len(self.runs)
        if run_count <= 1:
            errors = [math.nan] * len(self.histogram)
        else:
            jackknife = [[] for _ in self.histogram]
            for run in self.runs:
                histogram = [
                    total - contribution
                    for total, contribution in zip(self.histogram, run.histogram)
                ]
                leave_one_out, _ = self._normalized(
                    histogram,
                    self.weight_sum - run.weight_sum,
                    self.weighted_sigma_sum - run.weight_sum * run.sigma_gen,
                )
                for i, value in enumerate(leave_one_out):
                    jackknife[i].append(value)
            errors = []
            for estimates in jackknife:
                mean = sum(estimates) / run_count
                variance = (run_count - 1) / run_count * sum(
                    (value - mean) ** 2 for value in estimates
                )
                errors.append(math.sqrt(max(0.0, variance)))
        return SpectrumStats(
            values=values,
            standard_errors=errors,
            event_count=self.event_count,
            run_count=run_count,
            sigma_gen=sigma_gen,
            weight_sum=self.weight_sum,
        )


def parse_status_file(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in path.read_text(errors="replace").splitlines():
        if "=" in raw:
            key, value = raw.split("=", 1)
            out[key] = value
    return out


def successful_chunks(local_eos: Path, kind: str) -> list[int]:
    status_dir = local_eos / "status" / kind
    if not status_dir.exists():
        return []
    chunks = []
    for path in sorted(status_dir.glob("chunk_*.txt")):
        if parse_status_file(path).get("status") != "success":
            continue
        match = re.search(r"chunk_(\d+)\.txt$", path.name)
        if match:
            chunks.append(int(match.group(1)))
    return chunks


def copy_from_eos(eos_base: str, local_eos: Path) -> None:
    local_eos.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["xrdcp", "-r", "-f", f"root://eosuser.cern.ch/{eos_base.rstrip('/')}/", str(local_eos)],
        check=True,
    )


def eta_from_p(px: float, py: float, pz: float) -> float:
    pt = math.hypot(px, py)
    if pt == 0.0:
        return math.copysign(math.inf, pz)
    return math.asinh(pz / pt)


def charge_sign_from_label(label: int) -> float:
    return -1.0 if label in (2, 3) else 1.0


def parse_pythia_run(stream: BinaryIO, bins: list[float], eta_max: float) -> PythiaRun:
    widths = [bins[i + 1] - bins[i] for i in range(len(bins) - 1)]
    run_histogram = [0.0] * (len(bins) - 1)
    hist: list[float] | None = None
    weight: float | None = None
    current_cross: float | None = None
    final_cross: float | None = None
    weight_sum = 0.0
    event_count = 0

    def finish_event() -> None:
        nonlocal hist, weight, current_cross, final_cross, weight_sum, event_count
        if hist is None:
            return
        if weight is None or current_cross is None:
            raise ValueError("event is missing its PYTHIA weight or sigmaGen value")
        for i, value in enumerate(hist):
            run_histogram[i] += weight * value / widths[i]
        weight_sum += weight
        final_cross = current_cross
        event_count += 1
        hist = None
        weight = None
        current_cross = None

    for raw in stream:
        line = raw.decode("utf-8", errors="replace").strip()
        if not line:
            continue
        if line.startswith("# event"):
            finish_event()
            hist = [0.0] * (len(bins) - 1)
            weight = None
            current_cross = None
            continue
        if line.startswith("weight"):
            parts = line.split()
            weight = float(parts[1])
            if "cross" in parts:
                index = parts.index("cross")
                if index + 1 < len(parts):
                    current_cross = float(parts[index + 1])
            continue
        if line == "end":
            finish_event()
            continue
        parts = line.split()
        if hist is None or len(parts) != 6:
            continue
        px, py, pz = map(float, parts[:3])
        pdg = int(float(parts[4]))
        label = int(float(parts[5]))
        if label == -2 or abs(pdg) not in CHARGED_ABS:
            continue
        if abs(eta_from_p(px, py, pz)) > eta_max:
            continue
        pt = math.hypot(px, py)
        sign = charge_sign_from_label(label)
        for i in range(len(bins) - 1):
            if bins[i] <= pt < bins[i + 1]:
                hist[i] += sign
                break

    finish_event()
    if event_count == 0 or final_cross is None:
        raise ValueError("HYBRID output contains no complete PYTHIA events")
    if not math.isfinite(weight_sum) or weight_sum == 0.0:
        raise ValueError(f"invalid PYTHIA weightSum reconstructed from output: {weight_sum}")
    if not math.isfinite(final_cross):
        raise ValueError(f"invalid final PYTHIA sigmaGen value: {final_cross}")
    return PythiaRun(run_histogram, final_cross, weight_sum, event_count)


def read_summary(text: str) -> dict[str, str] | None:
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) < 2:
        return None
    return next(csv.DictReader(lines, delimiter="\t"))


def infer_variant(member_name: str, kind: str, summary: dict[str, str] | None) -> str:
    if kind == "pp":
        return "pp_reference"
    if summary is not None and summary.get("variant"):
        return summary["variant"]
    parent = PurePosixPath(member_name).parent.name
    return "with_prehydro" if parent.endswith("_prehydro") else "no_prehydro"


def scan_tar(
    tar_path: Path,
    kind: str,
    bins: list[float],
    eta_max: float,
) -> list[tuple[str, int | None, PythiaRun]]:
    found: list[tuple[str, int | None, PythiaRun]] = []
    with tarfile.open(tar_path, "r:gz") as tar:
        summaries: dict[str, dict[str, str]] = {}
        for member in tar.getmembers():
            if not member.name.endswith("/summary.tsv"):
                continue
            extracted = tar.extractfile(member)
            if extracted is None:
                continue
            summary = read_summary(extracted.read().decode("utf-8", errors="replace"))
            if summary is not None:
                summaries[str(PurePosixPath(member.name).parent)] = summary

        for member in tar.getmembers():
            if not member.name.endswith("/HYBRID_Hadrons.out"):
                continue
            extracted = tar.extractfile(member)
            if extracted is None:
                continue
            parent = str(PurePosixPath(member.name).parent)
            summary = summaries.get(parent)
            variant = infer_variant(member.name, kind, summary)
            hydro_index = None
            match = re.search(r"/hydro_(\d+)_", member.name)
            if match:
                hydro_index = int(match.group(1))
            found.append((variant, hydro_index, parse_pythia_run(extracted, bins, eta_max)))
    return found


def ratio_with_error(num: float, num_se: float, den: float, den_se: float) -> tuple[float, float]:
    if not math.isfinite(num) or not math.isfinite(den) or den == 0.0:
        return math.nan, math.nan
    ratio = num / den
    rel2 = 0.0
    if num != 0.0 and math.isfinite(num_se):
        rel2 += (num_se / num) ** 2
    if math.isfinite(den_se):
        rel2 += (den_se / den) ** 2
    return ratio, abs(ratio) * math.sqrt(rel2)


def write_variant_table(
    out_path: Path,
    bins: list[float],
    variant: str,
    aa: PythiaAggregate,
    pp: PythiaAggregate,
) -> None:
    aa_stats = aa.stats()
    pp_stats = pp.stats()
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "variant", "pt_low", "pt_high", "pt_center", "raa", "stat_err",
            "aa_events", "pp_events", "aa_spectrum", "aa_spectrum_stat_err",
            "pp_spectrum", "pp_spectrum_stat_err", "aa_runs", "pp_runs",
            "aa_sigma_gen", "pp_sigma_gen", "aa_weight_sum", "pp_weight_sum",
        ])
        for i in range(len(bins) - 1):
            raa, err = ratio_with_error(
                aa_stats.values[i],
                aa_stats.standard_errors[i],
                pp_stats.values[i],
                pp_stats.standard_errors[i],
            )
            writer.writerow([
                variant,
                f"{bins[i]:.8g}",
                f"{bins[i + 1]:.8g}",
                f"{0.5 * (bins[i] + bins[i + 1]):.8g}",
                f"{raa:.10g}",
                f"{err:.10g}",
                aa_stats.event_count,
                pp_stats.event_count,
                f"{aa_stats.values[i]:.10e}",
                f"{aa_stats.standard_errors[i]:.10e}",
                f"{pp_stats.values[i]:.10e}",
                f"{pp_stats.standard_errors[i]:.10e}",
                aa_stats.run_count,
                pp_stats.run_count,
                f"{aa_stats.sigma_gen:.10e}",
                f"{pp_stats.sigma_gen:.10e}",
                f"{aa_stats.weight_sum:.10e}",
                f"{pp_stats.weight_sum:.10e}",
            ])


def maybe_plot(out_dir: Path, overlay_tsv: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        (out_dir / "plot_error.txt").write_text(str(exc) + "\n")
        return

    rows = list(csv.DictReader(overlay_tsv.open(), delimiter="\t"))
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    styles = {
        "no_prehydro": {"label": "OO 0-5%, no prehydro", "color": "#0072B2", "marker": "o"},
        "with_prehydro": {
            "label": "OO 0-5%, public 2509.19430v2 Plan B",
            "color": "#D55E00",
            "marker": "s",
        },
    }
    for variant in AA_VARIANTS:
        selected = [row for row in rows if row["variant"] == variant]
        if not selected:
            continue
        x = [float(row["pt_center"]) for row in selected]
        y = [float(row["raa"]) for row in selected]
        ye = [float(row["stat_err"]) if row["stat_err"] != "nan" else math.nan for row in selected]
        style = styles[variant]
        ax.errorbar(x, y, yerr=ye, linewidth=1.8, linestyle="-", capsize=2.5, **style)
    ax.axhline(1.0, color="0.55", linewidth=1.0)
    ax.set_xlabel(r"charged hadron $p_T$ [GeV]")
    ax.set_ylabel(r"$R_{\rm AA}$")
    ax.set_title(r"O16+O16 5.36 TeV, 0--5%, charged $\pi/K/p$, $|\eta|<1$")
    ax.set_xlim(3.5, 112)
    ax.set_ylim(0.0, 1.5)
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "oo5360_c0_5_prehydro_overlay_raa.pdf")
    fig.savefig(out_dir / "oo5360_c0_5_prehydro_overlay_raa.png", dpi=180)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--local-eos", type=Path, required=True)
    parser.add_argument(
        "--pp-local-eos",
        type=Path,
        help="optional separate EOS snapshot containing the pp denominator",
    )
    parser.add_argument("--eos-base", help="Optional EOS path to xrdcp into --local-eos before analysis.")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bins", default=",".join(str(x) for x in DEFAULT_BINS))
    parser.add_argument("--eta-max", type=float, default=1.0)
    parser.add_argument(
        "--require-paired-aa",
        action="store_true",
        help="accept an AA chunk only when both variants have the required event count",
    )
    parser.add_argument(
        "--aa-events-per-chunk",
        type=int,
        default=0,
        help="required event count per AA variant when --require-paired-aa is set",
    )
    parser.add_argument("--copy", action="store_true")
    args = parser.parse_args()

    if args.copy:
        if not args.eos_base:
            raise ValueError("--copy requires --eos-base")
        copy_from_eos(args.eos_base, args.local_eos)

    bins = [float(x) for x in args.bins.split(",") if x]
    pp_local_eos = args.pp_local_eos or args.local_eos
    pp = PythiaAggregate(len(bins) - 1)
    aa = {variant: PythiaAggregate(len(bins) - 1) for variant in AA_VARIANTS}

    pp_missing_outputs = 0
    for chunk in successful_chunks(pp_local_eos, "pp"):
        tar_path = pp_local_eos / "outputs" / "pp" / f"chunk_{chunk}.tar.gz"
        if not tar_path.is_file():
            pp_missing_outputs += 1
            continue
        for _, _, run in scan_tar(tar_path, "pp", bins, args.eta_max):
            pp.add(run)

    aa_chunks_used = 0
    aa_chunks_skipped = 0
    aa_missing_outputs = 0
    for chunk in successful_chunks(args.local_eos, "aa"):
        tar_path = args.local_eos / "outputs" / "aa" / f"chunk_{chunk}.tar.gz"
        if not tar_path.is_file():
            aa_missing_outputs += 1
            continue
        chunk_runs = {variant: [] for variant in AA_VARIANTS}
        for variant, hydro_index, run in scan_tar(tar_path, "aa", bins, args.eta_max):
            if hydro_index not in (None, 0):
                continue
            if variant in aa:
                chunk_runs[variant].append(run)

        if args.require_paired_aa:
            expected = args.aa_events_per_chunk
            event_counts = [
                sum(run.event_count for run in chunk_runs[variant])
                for variant in AA_VARIANTS
            ]
            one_run_per_variant = all(
                len(chunk_runs[variant]) == 1 for variant in AA_VARIANTS
            )
            complete = (
                one_run_per_variant
                and (
                    all(count == expected for count in event_counts)
                    if expected > 0
                    else event_counts[0] > 0 and event_counts[0] == event_counts[1]
                )
            )
            if not complete:
                aa_chunks_skipped += 1
                continue

        aa_chunks_used += 1
        for variant in AA_VARIANTS:
            for run in chunk_runs[variant]:
                aa[variant].add(run)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    overlay_path = args.out_dir / "oo5360_c0_5_prehydro_overlay_raa.tsv"
    with overlay_path.open("w", newline="") as handle:
        writer = None
        for variant in AA_VARIANTS:
            tmp_path = args.out_dir / f"oo5360_c0_5_{variant}_raa.tsv"
            write_variant_table(tmp_path, bins, variant, aa[variant], pp)
            rows = list(csv.reader(tmp_path.open(), delimiter="\t"))
            if writer is None:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(rows[0])
            writer.writerows(rows[1:])

    maybe_plot(args.out_dir, overlay_path)
    print("normalization: PYTHIA aggregate sigmaGen/weightSum (PythiaParallel convention)")
    print(f"pp runs: {len(pp.runs)}")
    print(f"pp events: {pp.event_count}")
    print(f"pp status entries missing output archives: {pp_missing_outputs}")
    print(f"pp source: {pp_local_eos}")
    print(f"aa chunks used: {aa_chunks_used}")
    print(f"aa chunks skipped: {aa_chunks_skipped}")
    print(f"aa status entries missing output archives: {aa_missing_outputs}")
    for variant in AA_VARIANTS:
        print(f"{variant} aa runs: {len(aa[variant].runs)}")
        print(f"{variant} aa events: {aa[variant].event_count}")
    print(f"wrote {overlay_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
