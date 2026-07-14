#!/usr/bin/env python3
"""Build aligned no-prehydro/alpha=0.37 and no-prehydro/alpha=0.335 ROOT pairs."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import shlex
import subprocess
import sys

import convert_oo_paired_to_root as paired


def read_task_ids(path: Path) -> list[int]:
    values = [int(line) for line in path.read_text().splitlines() if line.strip()]
    if not values or values != sorted(set(values)):
        raise ValueError("accepted task IDs must be nonempty, sorted, and unique")
    return values


def status_is_success(path: Path, task_id: int) -> bool:
    if not path.is_file():
        return False
    values = paired.parse_key_value_text(path.read_text(errors="replace"))
    return (
        values.get("kind") == "aa"
        and values.get("task_id") == str(task_id)
        and values.get("status") == "success"
        and values.get("exit_code") == "0"
    )


def writer_command(
    writer: Path,
    output: Path,
    *,
    source_name: str,
    source_path: Path,
    raw_jet_pt_min: float,
    jet_abs_eta_max: float,
    z_cut: float,
    beta: float,
    match_dr_fraction: float,
) -> list[str]:
    return [
        str(writer),
        "--output",
        str(output),
        "--raw-jet-pt-min",
        str(raw_jet_pt_min),
        "--jet-abs-eta-max",
        str(jet_abs_eta_max),
        "--z-cut",
        str(z_cut),
        "--beta",
        str(beta),
        "--match-dr-fraction",
        str(match_dr_fraction),
        "--source",
        f"{source_name}\t{source_path}",
    ]


def convert(args: argparse.Namespace) -> dict[str, object]:
    reference = args.reference_source.resolve()
    alpha_source = args.alpha_source.resolve()
    manifest = args.manifest.resolve()
    accepted_path = args.accepted_task_ids.resolve()
    task_ids = read_task_ids(accepted_path)
    assignments = paired.load_aa_task_manifest(manifest)
    missing_manifest = set(task_ids) - set(assignments)
    if missing_manifest:
        raise ValueError("accepted task IDs are not contained in the task manifest")

    outputs = {
        "alpha037": args.output_alpha037.resolve(),
        "alpha0335": args.output_alpha0335.resolve(),
    }
    for output in outputs.values():
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists() and not args.overwrite:
            raise FileExistsError(f"refusing to overwrite {output}")
    audit_path = args.audit_output.resolve()
    summary_path = args.summary_output.resolve()
    for path in (audit_path, summary_path):
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.exists() and not args.overwrite:
            raise FileExistsError(f"refusing to overwrite {path}")

    writer = paired.build_writer(
        args.writer_source.resolve(),
        args.build_dir.resolve(),
        args.cxx,
        args.force_build,
    )
    partials = {
        key: output.with_name(f".{output.name}.partial")
        for key, output in outputs.items()
    }
    logs = {
        key: output.with_suffix(".writer.log") for key, output in outputs.items()
    }
    for path in [*partials.values(), *logs.values()]:
        path.unlink(missing_ok=True)

    log_streams = {key: path.open("w") for key, path in logs.items()}
    processes: dict[str, subprocess.Popen[bytes]] = {}
    started = datetime.now(timezone.utc)
    try:
        for key, label in (
            ("alpha037", "matched_no_vs_alpha037"),
            ("alpha0335", "matched_no_vs_alpha0335"),
        ):
            command = writer_command(
                writer,
                partials[key],
                source_name=label,
                source_path=reference if key == "alpha037" else alpha_source,
                raw_jet_pt_min=args.raw_jet_pt_min,
                jet_abs_eta_max=args.jet_abs_eta_max,
                z_cut=args.z_cut,
                beta=args.beta,
                match_dr_fraction=args.match_dr_fraction,
            )
            processes[key] = subprocess.Popen(
                command,
                stdin=subprocess.PIPE,
                stdout=log_streams[key],
                stderr=log_streams[key],
            )

        audit_fields = [
            "task_id",
            "pair_id",
            "seed",
            "hydro_index",
            "hydro_event_id",
            "hydro_ncoll",
            "hydro_payload_sha256",
            "event_weight",
            "sigma_gen",
            "hard_x",
            "hard_y",
            "hard_marker_count",
            "no_prehydro_particles",
            "alpha037_particles",
            "alpha0335_particles",
            "status",
        ]
        with audit_path.open("w", newline="") as audit_stream:
            audit = csv.DictWriter(
                audit_stream, fieldnames=audit_fields, delimiter="\t", lineterminator="\n"
            )
            audit.writeheader()
            for index, task_id in enumerate(task_ids, 1):
                reference_status = reference / "status/aa" / f"chunk_{task_id}.txt"
                alpha_status = alpha_source / "status/aa" / f"chunk_{task_id}.txt"
                if not status_is_success(reference_status, task_id):
                    raise paired.ArchiveValidationError(
                        f"task {task_id}: reference status is not strict success"
                    )
                if not status_is_success(alpha_status, task_id):
                    raise paired.ArchiveValidationError(
                        f"task {task_id}: alpha=0.335 status is not strict success"
                    )
                reference_archive = (
                    reference / "outputs/aa" / f"chunk_{task_id}.tar.gz"
                )
                alpha_archive = (
                    alpha_source / "outputs/aa" / f"chunk_{task_id}.tar.gz"
                )
                pair = paired.parse_paired_archive(
                    reference_archive, expected_chunk_id=task_id
                )
                alpha = paired.parse_prehydro_only_archive(
                    alpha_archive,
                    expected_chunk_id=task_id,
                    expected_alpha=args.expected_alpha,
                    expected_broadening_k=args.expected_broadening_k,
                )
                assignment = assignments[task_id]
                paired.validate_pair_task_assignment(
                    pair, task_id=task_id, assignment=assignment
                )
                paired.validate_prehydro_only_task_assignment(
                    alpha, task_id=task_id, assignment=assignment
                )
                paired.validate_hard_event_identity(
                    pair.no_prehydro, pair.with_prehydro
                )
                paired.validate_hard_event_identity(pair.no_prehydro, alpha.event)

                replacement = paired.PairedArchive(
                    seed=alpha.seed,
                    hydro_index=alpha.hydro_index,
                    hydro_event_id=alpha.hydro_event_id,
                    hydro_ncoll=alpha.hydro_ncoll,
                    hydro_payload_sha256=alpha.hydro_payload_sha256,
                    no_prehydro=pair.no_prehydro,
                    with_prehydro=alpha.event,
                )
                pair_id = paired.stable_pair_id(0, task_id)
                for key, payload in (("alpha037", pair), ("alpha0335", replacement)):
                    stream = processes[key].stdin
                    if stream is None:
                        raise RuntimeError(f"{key} ROOT writer has no input stream")
                    paired.write_pair(stream, pair_id, 0, task_id, payload)

                hard_markers = sum(
                    particle.raw_label == -2
                    for particle in pair.no_prehydro.particles
                )
                audit.writerow(
                    {
                        "task_id": task_id,
                        "pair_id": pair_id,
                        "seed": pair.seed,
                        "hydro_index": pair.hydro_index,
                        "hydro_event_id": pair.hydro_event_id,
                        "hydro_ncoll": pair.hydro_ncoll,
                        "hydro_payload_sha256": pair.hydro_payload_sha256,
                        "event_weight": f"{pair.no_prehydro.event_weight:.17g}",
                        "sigma_gen": f"{pair.no_prehydro.sigma_gen:.17g}",
                        "hard_x": f"{pair.no_prehydro.hard_x:.17g}",
                        "hard_y": f"{pair.no_prehydro.hard_y:.17g}",
                        "hard_marker_count": hard_markers,
                        "no_prehydro_particles": len(pair.no_prehydro.particles),
                        "alpha037_particles": len(pair.with_prehydro.particles),
                        "alpha0335_particles": len(alpha.event.particles),
                        "status": "PASS",
                    }
                )
                if index % args.progress_every == 0:
                    audit_stream.flush()
                    elapsed = (datetime.now(timezone.utc) - started).total_seconds()
                    print(
                        f"accepted_triplets={index}/{len(task_ids)} elapsed={elapsed:.1f}s",
                        file=sys.stderr,
                        flush=True,
                    )

        for process in processes.values():
            if process.stdin is not None:
                process.stdin.close()
        for key, process in processes.items():
            return_code = process.wait()
            if return_code != 0:
                raise RuntimeError(
                    f"{key} ROOT writer exited with {return_code}; inspect {logs[key]}"
                )
    except BaseException:
        for process in processes.values():
            if process.stdin is not None:
                try:
                    process.stdin.close()
                except (BrokenPipeError, OSError):
                    pass
            if process.poll() is None:
                process.terminate()
            process.wait()
        for path in partials.values():
            path.unlink(missing_ok=True)
        raise
    finally:
        for stream in log_streams.values():
            stream.close()

    root_totals = {
        key: paired.validate_root(path, len(task_ids))
        for key, path in partials.items()
    }
    for key, path in partials.items():
        path.replace(outputs[key])
    finished = datetime.now(timezone.utc)
    summary: dict[str, object] = {
        "status": "PASS",
        "schemaVersion": "oo-v3-two-aligned-pairs-v1",
        "startedUtc": started.isoformat(),
        "finishedUtc": finished.isoformat(),
        "elapsedSeconds": (finished - started).total_seconds(),
        "acceptedTriplets": len(task_ids),
        "expectedAlpha": args.expected_alpha,
        "expectedBroadeningK": args.expected_broadening_k,
        "identityChecks": [
            "task_id",
            "seed",
            "hydro_slot",
            "hydro_event_id",
            "hydro_ncoll",
            "hydro_payload_sha256",
            "event_number",
            "event_weight",
            "sigma_gen",
            "hard_x",
            "hard_y",
            "outgoing_hard_parton_markers",
        ],
        "taskIdList": str(accepted_path),
        "taskIdListSha256": paired.sha256(accepted_path),
        "taskManifest": str(manifest),
        "taskManifestSha256": paired.sha256(manifest),
        "referenceSource": str(reference),
        "alpha0335Source": str(alpha_source),
        "outputs": {
            key: {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": paired.sha256(path),
                "rootTotals": root_totals[key],
                "writerLog": str(logs[key]),
            }
            for key, path in outputs.items()
        },
        "writerBinary": str(writer),
        "rootVersion": paired.command_output(["root-config", "--version"]),
        "fastjetVersion": paired.command_output(["fastjet-config", "--version"]),
        "command": shlex.join(sys.argv),
    }
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-source", type=Path, required=True)
    parser.add_argument("--alpha-source", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--accepted-task-ids", type=Path, required=True)
    parser.add_argument("--output-alpha037", type=Path, required=True)
    parser.add_argument("--output-alpha0335", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument(
        "--writer-source",
        type=Path,
        default=Path(__file__).with_name("oo_root_tree_writer.cc"),
    )
    parser.add_argument("--cxx", default="c++")
    parser.add_argument("--force-build", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--progress-every", type=int, default=500)
    parser.add_argument("--expected-alpha", type=float, default=0.335)
    parser.add_argument("--expected-broadening-k", type=float, default=15.0)
    parser.add_argument("--raw-jet-pt-min", type=float, default=1.0)
    parser.add_argument("--jet-abs-eta-max", type=float, default=5.0)
    parser.add_argument("--z-cut", type=float, default=0.1)
    parser.add_argument("--beta", type=float, default=0.0)
    parser.add_argument("--match-dr-fraction", type=float, default=0.5)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.progress_every <= 0:
        raise ValueError("--progress-every must be positive")
    summary = convert(args)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
