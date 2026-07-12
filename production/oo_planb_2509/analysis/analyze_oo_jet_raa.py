#!/usr/bin/env python3
"""Calculate OO jet RAA from paired AA ROOT jets and streamed pp archives."""

from __future__ import annotations

import argparse
import awkward as ak
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tarfile

import numpy as np
import uproot


DEFAULT_BINS = [20, 25, 30, 40, 50, 65, 80, 110, 150, 220]
RADII = (0.1, 0.2, 0.4, 0.8)
RADIUS_DIGITS = {0.1: 1, 0.2: 2, 0.4: 4, 0.8: 8}
VARIANTS = ("noPrehydro", "withPrehydro")


def parse_status(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for line in path.read_text(errors="replace").splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            values[key] = value
    return values


def successful_pp_chunks(local_eos: Path) -> list[int]:
    chunks = []
    for path in sorted((local_eos / "status/pp").glob("chunk_*.txt")):
        if parse_status(path).get("status") != "success":
            continue
        chunks.append(int(path.stem.split("_")[-1]))
    return sorted(chunks)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def command_output(command: list[str]) -> str:
    return subprocess.run(command, check=True, text=True, capture_output=True).stdout.strip()


def build_helper(source: Path, build_dir: Path, cxx: str, force: bool = False) -> Path:
    build_dir.mkdir(parents=True, exist_ok=True)
    binary = build_dir / "oo_pp_jet_spectrum"
    if not force and binary.exists() and binary.stat().st_mtime_ns >= source.stat().st_mtime_ns:
        return binary
    flags = shlex.split(command_output(["fastjet-config", "--cxxflags"]))
    libraries = shlex.split(command_output(["fastjet-config", "--libs"]))
    temporary = binary.with_suffix(".tmp")
    subprocess.run(
        [
            cxx,
            "-O3",
            "-DNDEBUG",
            "-Wall",
            "-Wextra",
            "-Wpedantic",
            *flags,
            str(source),
            "-o",
            str(temporary),
            *libraries,
        ],
        check=True,
    )
    temporary.replace(binary)
    return binary


def pp_inventory(local_eos: Path) -> tuple[list[int], str, int]:
    chunks = successful_pp_chunks(local_eos)
    rows = []
    total_bytes = 0
    for chunk in chunks:
        archive = local_eos / "outputs/pp" / f"chunk_{chunk}.tar.gz"
        if not archive.is_file():
            raise FileNotFoundError(f"missing successful pp archive {archive}")
        size = archive.stat().st_size
        total_bytes += size
        rows.append((chunk, size))
    digest = hashlib.sha256(json.dumps(rows, separators=(",", ":")).encode()).hexdigest()
    return chunks, digest, total_bytes


def stream_pp_archives(
    local_eos: Path,
    chunks: list[int],
    helper: Path,
    output: Path,
    bins: np.ndarray,
    eta_max: float,
    log_path: Path,
) -> None:
    command = [
        str(helper),
        "--output",
        str(output),
        "--bins",
        ",".join(f"{value:g}" for value in bins),
        "--jet-abs-eta-max",
        f"{eta_max:g}",
        "--progress-every",
        "25",
    ]
    with log_path.open("wb") as log:
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=log, stderr=log)
        assert process.stdin is not None
        try:
            for chunk in chunks:
                archive = local_eos / "outputs/pp" / f"chunk_{chunk}.tar.gz"
                with tarfile.open(archive, "r:gz") as tar:
                    members = [
                        member
                        for member in tar.getmembers()
                        if member.isfile() and member.name.endswith("/HYBRID_Hadrons.out")
                    ]
                    if len(members) != 1:
                        raise ValueError(f"{archive}: expected one HYBRID_Hadrons.out, found {len(members)}")
                    extracted = tar.extractfile(members[0])
                    if extracted is None:
                        raise ValueError(f"{archive}: could not read {members[0].name}")
                    process.stdin.write(f"# OORUN {chunk}\n".encode())
                    shutil.copyfileobj(extracted, process.stdin, length=8 * 1024 * 1024)
                    process.stdin.write(b"\n# OOENDRUN\n")
        except Exception:
            process.stdin.close()
            process.terminate()
            process.wait()
            raise
        process.stdin.close()
        return_code = process.wait()
    if return_code != 0:
        raise RuntimeError(f"pp FastJet helper failed with {return_code}; inspect {log_path}")


def jackknife_error(leave_one_out: np.ndarray) -> np.ndarray:
    finite = np.all(np.isfinite(leave_one_out), axis=0)
    errors = np.full(leave_one_out.shape[1], np.nan, dtype=np.float64)
    if np.any(finite):
        values = leave_one_out[:, finite]
        means = np.mean(values, axis=0)
        count = values.shape[0]
        errors[finite] = np.sqrt(
            (count - 1.0) / count * np.sum((values - means) ** 2, axis=0)
        )
    return errors


def aggregate_spectrum(
    weighted_density_by_run: np.ndarray,
    weight_sum_by_run: np.ndarray,
    sigma_gen_by_run: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    if weighted_density_by_run.shape[0] != len(weight_sum_by_run):
        raise ValueError("spectrum rows and run metadata differ")
    weight_sum = float(np.sum(weight_sum_by_run))
    weighted_sigma_sum = float(np.sum(weight_sum_by_run * sigma_gen_by_run))
    if weight_sum <= 0.0 or weighted_sigma_sum <= 0.0:
        raise ValueError("invalid PYTHIA aggregate")
    totals = np.sum(weighted_density_by_run, axis=0)
    values = weighted_sigma_sum * totals / (weight_sum * weight_sum)
    if len(weight_sum_by_run) <= 1:
        return (
            values,
            np.full(weighted_density_by_run.shape[1], np.nan, dtype=np.float64),
            weight_sum,
            weighted_sigma_sum / weight_sum,
        )
    leave_weight_sum = weight_sum - weight_sum_by_run
    leave_weighted_sigma = weighted_sigma_sum - weight_sum_by_run * sigma_gen_by_run
    leave_totals = totals[np.newaxis, :] - weighted_density_by_run
    leave_values = (
        leave_weighted_sigma[:, np.newaxis]
        * leave_totals
        / leave_weight_sum[:, np.newaxis] ** 2
    )
    return values, jackknife_error(leave_values), weight_sum, weighted_sigma_sum / weight_sum


def independent_ratio(
    numerator: np.ndarray,
    numerator_error: np.ndarray,
    denominator: np.ndarray,
    denominator_error: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    ratio = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan),
        where=denominator != 0.0,
    )
    relative_variance = np.zeros_like(ratio)
    np.divide(
        numerator_error**2,
        numerator**2,
        out=relative_variance,
        where=numerator != 0.0,
    )
    denominator_term = np.zeros_like(ratio)
    np.divide(
        denominator_error**2,
        denominator**2,
        out=denominator_term,
        where=denominator != 0.0,
    )
    return ratio, np.abs(ratio) * np.sqrt(relative_variance + denominator_term)


def load_pp_cache(
    path: Path, bins: np.ndarray
) -> tuple[
    dict[float, np.ndarray],
    dict[float, np.ndarray],
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    run_ids = sorted({int(row["run_id"]) for row in rows})
    run_index = {run_id: index for index, run_id in enumerate(run_ids)}
    matrices = {
        radius: np.full((len(run_ids), len(bins) - 1), np.nan, dtype=np.float64)
        for radius in RADII
    }
    raw_matrices = {
        radius: np.zeros((len(run_ids), len(bins) - 1), dtype=np.int64)
        for radius in RADII
    }
    event_counts = np.zeros(len(run_ids), dtype=np.int64)
    weight_sums = np.full(len(run_ids), np.nan, dtype=np.float64)
    sigma_gen = np.full(len(run_ids), np.nan, dtype=np.float64)
    for row in rows:
        radius = float(row["radius"])
        if radius not in matrices:
            raise ValueError(f"unexpected pp radius {radius}")
        index = run_index[int(row["run_id"])]
        bin_index = int(row["bin_index"])
        if not math.isclose(float(row["pt_low"]), bins[bin_index]) or not math.isclose(
            float(row["pt_high"]), bins[bin_index + 1]
        ):
            raise ValueError("pp cache binning does not match request")
        metadata = (int(row["event_count"]), float(row["weight_sum"]), float(row["sigma_gen"]))
        if event_counts[index] == 0:
            event_counts[index], weight_sums[index], sigma_gen[index] = metadata
        elif (
            event_counts[index] != metadata[0]
            or not math.isclose(weight_sums[index], metadata[1])
            or not math.isclose(sigma_gen[index], metadata[2])
        ):
            raise ValueError("inconsistent pp metadata within one run")
        if np.isfinite(matrices[radius][index, bin_index]):
            raise ValueError("duplicate pp cache row")
        matrices[radius][index, bin_index] = float(row["weighted_density"])
        raw_matrices[radius][index, bin_index] = int(row["raw_jet_count"])
    for radius in RADII:
        if not np.all(np.isfinite(matrices[radius])):
            raise ValueError(f"incomplete pp cache for R={radius}")
    if not np.all(event_counts > 0) or not np.all(weight_sums > 0.0) or not np.all(sigma_gen > 0.0):
        raise ValueError("invalid pp run metadata")
    return matrices, raw_matrices, weight_sums, sigma_gen, event_counts


def aa_histogram_matrix(
    pt: ak.Array,
    eta: ak.Array,
    weights: np.ndarray,
    bins: np.ndarray,
    eta_max: float,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    selected = np.isfinite(pt) & np.isfinite(eta) & (np.abs(eta) < eta_max)
    values = pt[selected]
    lengths = ak.to_numpy(ak.num(values, axis=1)).astype(np.int64)
    flat = ak.to_numpy(ak.flatten(values)).astype(np.float64)
    event_indices = np.repeat(np.arange(len(weights), dtype=np.int64), lengths)
    indices = np.searchsorted(bins, flat, side="left") - 1
    accepted = (indices >= 0) & (indices < len(bins) - 1)
    raw = np.zeros((len(weights), len(bins) - 1), dtype=np.float64)
    np.add.at(raw, (event_indices[accepted], indices[accepted]), 1.0)
    weighted_density = raw * weights[:, np.newaxis] / np.diff(bins)[np.newaxis, :]
    return (
        weighted_density,
        raw,
        int(np.count_nonzero(indices < 0)),
        int(np.count_nonzero(indices >= len(bins) - 1)),
    )


def plot_results(rows: list[dict[str, object]], out_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    styles = {
        "noPrehydro": {"label": "no prehydro", "color": "#0072B2", "marker": "o"},
        "withPrehydro": {"label": "Plan B prehydro", "color": "#D55E00", "marker": "s"},
    }

    def draw(ax, radius: float) -> None:
        for variant in VARIANTS:
            selected = [row for row in rows if row["variant"] == variant and row["radius"] == radius]
            x = np.array([row["pt_center"] for row in selected], dtype=float)
            low = np.array([row["pt_low"] for row in selected], dtype=float)
            high = np.array([row["pt_high"] for row in selected], dtype=float)
            ax.errorbar(
                x,
                np.array([row["raa"] for row in selected], dtype=float),
                xerr=[x - low, high - x],
                yerr=np.array([row["stat_error"] for row in selected], dtype=float),
                linewidth=1.5,
                capsize=2.2,
                **styles[variant],
            )
        ax.axhline(1.0, color="0.5", linewidth=1.0)
        ax.set_xscale("log")
        ax.set_xlim(19, 230)
        ax.set_ylim(0.0, 1.4)
        ax.set_title(rf"anti-$k_T$ $R={radius:.1f}$")
        ax.set_xlabel(r"jet $p_T$ [GeV]")
        ax.grid(alpha=0.22)

    outputs = []
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.4), sharey=True)
    for ax, radius in zip(axes.flat, RADII):
        draw(ax, radius)
    axes[0, 0].set_ylabel(r"jet $R_{\rm AA}$")
    axes[1, 0].set_ylabel(r"jet $R_{\rm AA}$")
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(r"O16+O16 5.36 TeV, 0--5%, $|\eta_{\rm jet}|<2$ (provisional)")
    fig.tight_layout()
    for suffix in ("pdf", "png"):
        path = out_dir / f"oo5360_c0_5_jet_raa_R01020408.{suffix}"
        fig.savefig(path, dpi=180 if suffix == "png" else None)
        outputs.append(path)
    plt.close(fig)

    for radius in RADII:
        fig, ax = plt.subplots(figsize=(6.6, 4.8))
        draw(ax, radius)
        ax.set_ylabel(r"jet $R_{\rm AA}$")
        ax.legend(frameon=False, fontsize=9)
        fig.suptitle(r"O16+O16 5.36 TeV, 0--5%, $|\eta_{\rm jet}|<2$ (provisional)")
        fig.tight_layout()
        for suffix in ("pdf", "png"):
            path = out_dir / f"oo5360_c0_5_jet_raa_R{RADIUS_DIGITS[radius]:02d}.{suffix}"
            fig.savefig(path, dpi=180 if suffix == "png" else None)
            outputs.append(path)
        plt.close(fig)
    return outputs


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--pp-local-eos", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--pp-cache", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--bins", default=",".join(str(value) for value in DEFAULT_BINS))
    parser.add_argument("--jet-abs-eta-max", type=float, default=2.0)
    parser.add_argument("--expected-aa-events", type=int)
    parser.add_argument("--expected-pp-runs", type=int, default=1000)
    parser.add_argument("--expected-pp-events", type=int, default=1_000_000)
    parser.add_argument("--helper-source", type=Path, default=Path(__file__).with_name("oo_pp_jet_spectrum.cc"))
    parser.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--force-pp", action="store_true")
    args = parser.parse_args()

    bins = np.array([float(value) for value in args.bins.split(",") if value], dtype=np.float64)
    if len(bins) < 2 or not np.all(np.diff(bins) > 0.0):
        parser.error("--bins must be strictly increasing")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.pp_cache.parent.mkdir(parents=True, exist_ok=True)

    chunks, inventory_sha, archive_bytes = pp_inventory(args.pp_local_eos)
    helper_source_sha = sha256(args.helper_source)
    fastjet_version = command_output(["fastjet-config", "--version"])
    cache_metadata_path = args.pp_cache.with_suffix(args.pp_cache.suffix + ".json")
    requested_cache = {
        "ppLocalEos": str(args.pp_local_eos.resolve()),
        "ppInventorySha256": inventory_sha,
        "ppRuns": len(chunks),
        "binsGeV": bins.tolist(),
        "jetAbsEtaMax": args.jet_abs_eta_max,
        "radii": list(RADII),
        "helperSourceSha256": helper_source_sha,
        "fastjetVersion": fastjet_version,
    }
    use_cache = False
    if args.pp_cache.is_file() and cache_metadata_path.is_file() and not args.force_pp:
        cache_metadata = json.loads(cache_metadata_path.read_text())
        use_cache = all(cache_metadata.get(key) == value for key, value in requested_cache.items())
    helper = build_helper(args.helper_source, args.build_dir, args.cxx)
    if not use_cache:
        temporary = args.pp_cache.with_suffix(args.pp_cache.suffix + ".tmp")
        log_path = args.out_dir / "pp_fastjet.log"
        stream_pp_archives(
            args.pp_local_eos, chunks, helper, temporary, bins, args.jet_abs_eta_max, log_path
        )
        temporary.replace(args.pp_cache)
        cache_metadata_path.write_text(json.dumps(requested_cache, indent=2, sort_keys=True) + "\n")

    (
        pp_matrices,
        pp_raw_matrices,
        pp_weight_sums,
        pp_sigma_gen,
        pp_event_counts,
    ) = load_pp_cache(args.pp_cache, bins)
    if len(pp_weight_sums) != args.expected_pp_runs:
        raise RuntimeError(f"accepted {len(pp_weight_sums)} pp runs, expected {args.expected_pp_runs}")
    if int(np.sum(pp_event_counts)) != args.expected_pp_events:
        raise RuntimeError(
            f"accepted {int(np.sum(pp_event_counts))} pp events, expected {args.expected_pp_events}"
        )

    with uproot.open(args.input_root) as root_file:
        pairs = root_file["Pairs"].arrays(["pairId", "eventWeight", "sigmaGen", "seed"], library="np")
        pair_ids = pairs["pairId"]
        weights = pairs["eventWeight"].astype(np.float64)
        sigma_gen = pairs["sigmaGen"].astype(np.float64)
        if args.expected_aa_events is not None and len(pair_ids) != args.expected_aa_events:
            raise RuntimeError(f"accepted {len(pair_ids)} AA pairs, expected {args.expected_aa_events}")
        if len(np.unique(pair_ids)) != len(pair_ids) or len(np.unique(pairs["seed"])) != len(pair_ids):
            raise ValueError("AA ROOT contains duplicate pair IDs or seeds")
        aa_matrices: dict[str, dict[float, np.ndarray]] = {variant: {} for variant in VARIANTS}
        aa_raw: dict[str, dict[float, np.ndarray]] = {variant: {} for variant in VARIANTS}
        aa_flow: dict[str, dict[float, tuple[int, int]]] = {variant: {} for variant in VARIANTS}
        for variant in VARIANTS:
            branches = ["pairId", "eventWeight", "sigmaGen"]
            for radius in RADII:
                digit = RADIUS_DIGITS[radius]
                branches.extend([f"jet{digit}Pt", f"jet{digit}Eta"])
            arrays = root_file[f"{variant}/Jets"].arrays(branches, library="ak")
            if not np.array_equal(ak.to_numpy(arrays.pairId), pair_ids):
                raise ValueError(f"{variant} pair ordering differs from Pairs")
            for radius in RADII:
                digit = RADIUS_DIGITS[radius]
                matrix, raw, underflow, overflow = aa_histogram_matrix(
                    arrays[f"jet{digit}Pt"], arrays[f"jet{digit}Eta"], weights, bins, args.jet_abs_eta_max
                )
                aa_matrices[variant][radius] = matrix
                aa_raw[variant][radius] = raw
                aa_flow[variant][radius] = (underflow, overflow)

    pp_spectra = {}
    for radius in RADII:
        pp_spectra[radius] = aggregate_spectrum(
            pp_matrices[radius], pp_weight_sums, pp_sigma_gen
        )

    output_rows: list[dict[str, object]] = []
    for radius in RADII:
        pp_values, pp_errors, pp_weight_sum, pp_sigma_merged = pp_spectra[radius]
        for variant in VARIANTS:
            aa_values, aa_errors, aa_weight_sum, aa_sigma_merged = aggregate_spectrum(
                aa_matrices[variant][radius], weights, sigma_gen
            )
            raa, raa_errors = independent_ratio(aa_values, aa_errors, pp_values, pp_errors)
            for index in range(len(bins) - 1):
                output_rows.append(
                    {
                        "variant": variant,
                        "radius": radius,
                        "pt_low": bins[index],
                        "pt_high": bins[index + 1],
                        "pt_center": math.sqrt(bins[index] * bins[index + 1]),
                        "raa": raa[index],
                        "stat_error": raa_errors[index],
                        "aa_spectrum_mb_per_gev": aa_values[index],
                        "aa_stat_error_mb_per_gev": aa_errors[index],
                        "pp_spectrum_mb_per_gev": pp_values[index],
                        "pp_stat_error_mb_per_gev": pp_errors[index],
                        "aa_raw_jets": int(np.sum(aa_raw[variant][radius][:, index])),
                        "pp_raw_jets": int(np.sum(pp_raw_matrices[radius][:, index])),
                        "aa_events": len(pair_ids),
                        "pp_events": int(np.sum(pp_event_counts)),
                        "aa_weight_sum": aa_weight_sum,
                        "pp_weight_sum": pp_weight_sum,
                        "aa_sigma_merged_mb": aa_sigma_merged,
                        "pp_sigma_merged_mb": pp_sigma_merged,
                    }
                )

    table_path = args.out_dir / "oo5360_c0_5_jet_raa_R01020408.tsv"
    fields = list(output_rows[0])
    with table_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in output_rows:
            writer.writerow(
                {
                    key: f"{value:.12e}" if isinstance(value, float) else value
                    for key, value in row.items()
                }
            )

    plots = plot_results(output_rows, args.out_dir)
    metadata = {
        "status": "PASS",
        "publicationStatus": "PROVISIONAL_DIAGNOSTIC_NOT_UNBIASED",
        "selection": f"anti-kt E-scheme; corrected jet pT; abs(jet eta)<{args.jet_abs_eta_max:g}",
        "radii": list(RADII),
        "binsGeV": bins.tolist(),
        "aaInputRoot": str(args.input_root.resolve()),
        "aaInputRootBytes": args.input_root.stat().st_size,
        "aaInputRootSha256": sha256(args.input_root),
        "aaEvents": len(pair_ids),
        "ppLocalEos": str(args.pp_local_eos.resolve()),
        "ppRuns": len(pp_weight_sums),
        "ppEvents": int(np.sum(pp_event_counts)),
        "ppArchiveBytes": archive_bytes,
        "ppInventorySha256": inventory_sha,
        "ppCache": str(args.pp_cache.resolve()),
        "ppCacheSha256": sha256(args.pp_cache),
        "ppRawJetsAllRadiiInRange": int(
            sum(np.sum(values) for values in pp_raw_matrices.values())
        ),
        "normalization": "PythiaParallel merged sigmaGen/sum(weight), applied once; no second Ncoll factor",
        "uncertainty": "delete-one-AA-event and delete-one-pp-run jackknives, combined independently",
        "helperSource": str(args.helper_source.resolve()),
        "helperSourceSha256": helper_source_sha,
        "fastjetVersion": fastjet_version,
        "sourceRepositoryHead": command_output(
            ["git", "-C", str(Path(__file__).resolve().parents[3]), "rev-parse", "HEAD"]
        ),
        "table": str(table_path),
        "plots": [str(path) for path in plots],
        "aaFlowCounts": {
            variant: {
                f"R{radius:.1f}": {"belowRange": aa_flow[variant][radius][0], "aboveRange": aa_flow[variant][radius][1]}
                for radius in RADII
            }
            for variant in VARIANTS
        },
    }
    metadata_path = args.out_dir / "oo5360_c0_5_jet_raa_R01020408_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    marker = args.out_dir / "jet_raa_complete.txt"
    marker.write_text(
        f"status=PASS\naa_events={len(pair_ids)}\npp_events={int(np.sum(pp_event_counts))}\n"
        f"table={table_path}\nmetadata={metadata_path}\n"
    )
    print(marker.read_text(), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
