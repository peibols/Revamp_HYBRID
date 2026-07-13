#!/usr/bin/env python3
"""Monitor and repair a task-matched, prehydro-only OO campaign."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import math
from pathlib import Path
import subprocess
import tarfile
import time

from supervise_oo_v2 import (
    log,
    parse_status,
    parse_tsv_member,
    query_jobs,
    remove_terminal_holds,
    submit_retry,
    sync_eos,
)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--campaign", required=True)
    parser.add_argument("--task-manifest", type=Path, required=True)
    parser.add_argument("--target", type=int, default=100_000)
    parser.add_argument("--expected-alpha", type=float, required=True)
    parser.add_argument("--expected-broadening-k", type=float, default=15.0)
    parser.add_argument("--sleep-s", type=int, default=900)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--cernctl", default="/data/yjlee/cernLxplus/cernctl")
    parser.add_argument("--cern-remote", default="lxplus")
    parser.add_argument("--schedd", default="bigbird101.cern.ch")
    parser.add_argument("--afs-work", type=Path, required=True)
    parser.add_argument("--max-jobs-per-submit", type=int, default=10_000)
    parser.add_argument("--retry-timeout-s", type=int, default=72_000)
    args = parser.parse_args()
    if args.target <= 0 or args.sleep_s <= 0 or args.max_jobs_per_submit <= 0:
        parser.error("--target, --sleep-s, and --max-jobs-per-submit must be positive")
    if args.retry_timeout_s <= 0:
        parser.error("--retry-timeout-s must be positive")
    for name in ("expected_alpha", "expected_broadening_k"):
        value = getattr(args, name)
        if not math.isfinite(value) or value < 0:
            parser.error(f"--{name.replace('_', '-')} must be finite and nonnegative")
    return args


def load_assignments(path: Path, target: int) -> dict[int, dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != target:
        raise ValueError(f"{path}: expected {target} tasks, found {len(rows)}")
    assignments = {int(row["task_id"]): row for row in rows}
    if set(assignments) != set(range(target)):
        raise ValueError(f"{path}: task IDs must be exactly 0..{target - 1}")
    return assignments


def float_matches(raw: str | None, expected: float) -> bool:
    try:
        value = float(raw) if raw is not None else math.nan
    except ValueError:
        return False
    return math.isclose(value, expected, rel_tol=0.0, abs_tol=1e-12)


def row_matches_assignment(
    row: dict[str, str],
    assignment: dict[str, str],
    *,
    expected_alpha: float,
    expected_broadening_k: float,
) -> bool:
    expected = {
        "task_id": assignment["task_id"],
        "seed": assignment["hard_seed"],
        "centrality": "C0-5",
        "hydro_event_id": assignment["hydro_event_id"],
        "hydro_ncoll": assignment["hydro_ncoll"],
        "hydro_payload_sha256": assignment["hydro_payload_sha256"],
        "variant": "with_prehydro",
        "use_prehydro": "1",
        "returncode": "0",
        "timeout": "0",
    }
    if any(row.get(key) != value for key, value in expected.items()):
        return False
    hydro_slot = row.get("hydro_slot", row.get("hydro_index"))
    if hydro_slot != assignment["hydro_slot"]:
        return False
    return float_matches(
        row.get("energy_loss_alpha"), expected_alpha
    ) and float_matches(row.get("broadening_k"), expected_broadening_k)


def status_matches_assignment(
    status: dict[str, str],
    assignment: dict[str, str],
    *,
    expected_alpha: float,
    expected_broadening_k: float,
) -> bool:
    expected = {
        "kind": "aa",
        "task_id": assignment["task_id"],
        "status": "success",
        "exit_code": "0",
        "run_prehydro_pair": "false",
        "run_prehydro_only": "true",
        "hydro_slot": assignment["hydro_slot"],
        "hydro_event_id": assignment["hydro_event_id"],
        "hydro_ncoll": assignment["hydro_ncoll"],
        "hydro_dir": assignment["hydro_dir"],
        "hydro_payload_key": assignment["hydro_payload_key"],
        "hydro_payload_sha256": assignment["hydro_payload_sha256"],
    }
    if any(status.get(key) != value for key, value in expected.items()):
        return False
    return float_matches(
        status.get("prehydro_alpha"), expected_alpha
    ) and float_matches(status.get("broadening_k"), expected_broadening_k)


def prehydro_only_archive_is_complete(
    path: Path,
    assignment: dict[str, str],
    *,
    expected_alpha: float,
    expected_broadening_k: float,
) -> bool:
    task_id = int(assignment["task_id"])
    task = f"task_{task_id:05d}"
    required_suffixes = {
        "summary": f"/{task}_prehydro/summary.tsv",
        "hadrons": f"/{task}_prehydro/HYBRID_Hadrons.out",
        "input": f"/{task}_prehydro/hybrid_input.dat",
        "variant_summary": f"/{task}_prehydro_only_summary.tsv",
    }
    forbidden_suffixes = (
        f"/{task}/summary.tsv",
        f"/{task}/HYBRID_Hadrons.out",
        f"/{task}_pair_summary.tsv",
    )
    try:
        with tarfile.open(path, "r:gz") as archive:
            members = archive.getmembers()
            if any(
                member.name.endswith(suffix)
                for member in members
                for suffix in forbidden_suffixes
            ):
                return False
            matched: dict[str, tarfile.TarInfo] = {}
            for key, suffix in required_suffixes.items():
                candidates = [member for member in members if member.name.endswith(suffix)]
                if len(candidates) != 1 or not candidates[0].isfile():
                    return False
                matched[key] = candidates[0]
            if matched["hadrons"].size <= 0 or matched["input"].size <= 0:
                return False
            for key in ("summary", "variant_summary"):
                rows = parse_tsv_member(archive, matched[key].name)
                if len(rows) != 1 or not row_matches_assignment(
                    rows[0],
                    assignment,
                    expected_alpha=expected_alpha,
                    expected_broadening_k=expected_broadening_k,
                ):
                    return False
    except (KeyError, OSError, tarfile.TarError, UnicodeError, ValueError):
        return False
    return True


def ready_ids(
    local_eos: Path,
    assignments: dict[int, dict[str, str]],
    *,
    expected_alpha: float,
    expected_broadening_k: float,
) -> set[int]:
    ready: set[int] = set()
    for task_id, assignment in assignments.items():
        status_path = local_eos / f"status/aa/chunk_{task_id}.txt"
        output_path = local_eos / f"outputs/aa/chunk_{task_id}.tar.gz"
        if not status_path.is_file() or not output_path.is_file():
            continue
        if not status_matches_assignment(
            parse_status(status_path),
            assignment,
            expected_alpha=expected_alpha,
            expected_broadening_k=expected_broadening_k,
        ):
            continue
        if prehydro_only_archive_is_complete(
            output_path,
            assignment,
            expected_alpha=expected_alpha,
            expected_broadening_k=expected_broadening_k,
        ):
            ready.add(task_id)
    return ready


def supervise_once(args: argparse.Namespace) -> bool:
    assignments = load_assignments(args.task_manifest, args.target)
    local_eos = args.work / "eos_snapshot"
    seen_path = args.work / "prehydro_only_queue_seen.txt"
    states, active_aa, held_aa = query_jobs(args)
    if active_aa and not seen_path.exists():
        seen_path.write_text(
            f"date={datetime.now().astimezone().isoformat(timespec='seconds')}\n"
        )
    sync_eos(args, local_eos)
    ready = ready_ids(
        local_eos,
        assignments,
        expected_alpha=args.expected_alpha,
        expected_broadening_k=args.expected_broadening_k,
    )
    log(
        f"Condor states={dict(sorted(states.items()))}; "
        f"active_AA={len(active_aa)} strict_ready={len(ready)}/{args.target}"
    )
    active_aa -= remove_terminal_holds(args, held_aa, ready)
    if active_aa or not seen_path.exists():
        return False

    missing = set(assignments) - ready
    if missing:
        submit_retry(args, missing)
        return False

    marker = args.work / "prehydro_only_100k_strict_complete.txt"
    marker.write_text(
        f"date={datetime.now().astimezone().isoformat(timespec='seconds')}\n"
        f"accepted_variants={len(ready)}\n"
        f"prehydro_alpha={args.expected_alpha}\n"
        f"task_manifest={args.task_manifest}\n"
    )
    log("prehydro-only strict target complete")
    return True


def main() -> int:
    args = parse_args()
    while True:
        try:
            if supervise_once(args):
                return 0
        except (OSError, RuntimeError, ValueError, subprocess.SubprocessError) as error:
            log(f"supervisor iteration failed: {type(error).__name__}: {error}")
        if args.once:
            return 0
        time.sleep(args.sleep_s)


if __name__ == "__main__":
    raise SystemExit(main())
