#!/usr/bin/env python3
"""Analyze paired OO no-prehydro/prehydro AA outputs against a common pp sample."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path, PurePosixPath
import re
import subprocess
import tarfile
from typing import BinaryIO, Iterable


DEFAULT_BINS = [4, 5, 6, 8, 10, 15, 20, 30, 40, 60, 80, 110]
CHARGED_ABS = {211, 321, 2212}
AA_VARIANTS = ["no_prehydro", "with_prehydro"]


class RunningSpectrum:
    def __init__(self, n_bins: int) -> None:
        self.n = 0
        self.mean = [0.0] * n_bins
        self.m2 = [0.0] * n_bins

    def add(self, values: list[float]) -> None:
        self.n += 1
        for i, value in enumerate(values):
            delta = value - self.mean[i]
            self.mean[i] += delta / self.n
            self.m2[i] += delta * (value - self.mean[i])

    def mean_and_se(self) -> tuple[list[float], list[float], int]:
        if self.n == 0:
            return [math.nan] * len(self.mean), [math.nan] * len(self.mean), 0
        if self.n == 1:
            return list(self.mean), [math.nan] * len(self.mean), self.n
        return list(self.mean), [math.sqrt(x / (self.n - 1) / self.n) for x in self.m2], self.n


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


def parse_event_spectra(stream: BinaryIO, bins: list[float], eta_max: float) -> Iterable[list[float]]:
    widths = [bins[i + 1] - bins[i] for i in range(len(bins) - 1)]
    events: list[list[float]] = []
    cross_values: list[float] = []
    hist: list[float] | None = None
    weight = 1.0
    current_cross: float | None = None

    def finish_event() -> None:
        nonlocal hist, current_cross
        if hist is None:
            return
        events.append([value / widths[i] for i, value in enumerate(hist)])
        if current_cross is not None:
            cross_values.append(current_cross)
        hist = None
        current_cross = None

    for raw in stream:
        line = raw.decode("utf-8", errors="replace").strip()
        if not line:
            continue
        if line.startswith("# event"):
            finish_event()
            hist = [0.0] * (len(bins) - 1)
            weight = 1.0
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
                hist[i] += sign * weight
                break

    finish_event()
    average_cross = sum(cross_values) / len(cross_values) if cross_values else 1.0
    for event in events:
        yield [average_cross * value for value in event]


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
) -> list[tuple[str, int | None, list[float]]]:
    found: list[tuple[str, int | None, list[float]]] = []
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
            for event in parse_event_spectra(extracted, bins, eta_max):
                found.append((variant, hydro_index, event))
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
    aa: RunningSpectrum,
    pp: RunningSpectrum,
) -> None:
    aa_mean, aa_se, aa_events = aa.mean_and_se()
    pp_mean, pp_se, pp_events = pp.mean_and_se()
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "variant", "pt_low", "pt_high", "pt_center", "raa", "stat_err",
            "aa_events", "pp_events", "aa_spectrum", "aa_spectrum_stat_err",
            "pp_spectrum", "pp_spectrum_stat_err",
        ])
        for i in range(len(bins) - 1):
            raa, err = ratio_with_error(aa_mean[i], aa_se[i], pp_mean[i], pp_se[i])
            writer.writerow([
                variant,
                f"{bins[i]:.8g}",
                f"{bins[i + 1]:.8g}",
                f"{0.5 * (bins[i] + bins[i + 1]):.8g}",
                f"{raa:.10g}",
                f"{err:.10g}",
                aa_events,
                pp_events,
                f"{aa_mean[i]:.10e}",
                f"{aa_se[i]:.10e}",
                f"{pp_mean[i]:.10e}",
                f"{pp_se[i]:.10e}",
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
    pp = RunningSpectrum(len(bins) - 1)
    aa = {variant: RunningSpectrum(len(bins) - 1) for variant in AA_VARIANTS}

    for chunk in successful_chunks(pp_local_eos, "pp"):
        tar_path = pp_local_eos / "outputs" / "pp" / f"chunk_{chunk}.tar.gz"
        for _, _, event in scan_tar(tar_path, "pp", bins, args.eta_max):
            pp.add(event)

    aa_chunks_used = 0
    aa_chunks_skipped = 0
    for chunk in successful_chunks(args.local_eos, "aa"):
        tar_path = args.local_eos / "outputs" / "aa" / f"chunk_{chunk}.tar.gz"
        chunk_events = {variant: [] for variant in AA_VARIANTS}
        for variant, hydro_index, event in scan_tar(tar_path, "aa", bins, args.eta_max):
            if hydro_index not in (None, 0):
                continue
            if variant in aa:
                chunk_events[variant].append(event)

        if args.require_paired_aa:
            counts = [len(chunk_events[variant]) for variant in AA_VARIANTS]
            expected = args.aa_events_per_chunk
            complete = all(count == expected for count in counts) if expected > 0 else counts[0] > 0 and counts[0] == counts[1]
            if not complete:
                aa_chunks_skipped += 1
                continue

        aa_chunks_used += 1
        for variant in AA_VARIANTS:
            for event in chunk_events[variant]:
                aa[variant].add(event)

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
    print(f"pp events: {pp.n}")
    print(f"pp source: {pp_local_eos}")
    print(f"aa chunks used: {aa_chunks_used}")
    print(f"aa chunks skipped: {aa_chunks_skipped}")
    for variant in AA_VARIANTS:
        print(f"{variant} aa events: {aa[variant].n}")
    print(f"wrote {overlay_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
