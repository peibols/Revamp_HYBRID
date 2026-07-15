#!/usr/bin/env python3
"""Monitor and strictly validate the paired Moliere-on OO campaign."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import math
from pathlib import Path
import re
import subprocess
import tarfile
import time

from supervise_oo_v2 import (
    campaign_constraint,
    log,
    parse_status,
    parse_tsv_member,
    query_jobs,
    remove_terminal_holds,
    submit_retry,
    sync_eos,
)


def parse_args() -> argparse.Namespace:
    root = Path("/raid5/data/yjlee/hybrid_dev")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--work",
        type=Path,
        default=root / "test/oo5360_v4_moliere_pair_alpha0335_100k_20260715",
    )
    parser.add_argument(
        "--campaign",
        default=(
            "hybrid_oo5360_c0_5_500hydro_moliere_paired_"
            "alpha0335_v4_100kAA_20260715"
        ),
    )
    parser.add_argument("--task-manifest", type=Path, required=True)
    parser.add_argument("--task-id-start", type=int, default=0)
    parser.add_argument("--target", type=int, default=100_000)
    parser.add_argument("--expected-alpha", type=float, default=0.335)
    parser.add_argument("--expected-broadening-k", type=float, default=15.0)
    parser.add_argument("--expected-tables-sha256", required=True)
    parser.add_argument("--job-priority", type=int, default=50)
    parser.add_argument("--sleep-s", type=int, default=900)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--cernctl", default="/data/yjlee/cernLxplus/cernctl")
    parser.add_argument("--cern-remote", default="lxplus")
    parser.add_argument("--schedd", default="bigbird101.cern.ch")
    parser.add_argument(
        "--afs-work",
        type=Path,
        default=Path("/afs/cern.ch/user/y/yjlee/oo_v4_moliere_pair_alpha0335_100k_20260715"),
    )
    parser.add_argument("--max-jobs-per-submit", type=int, default=10_000)
    parser.add_argument("--retry-timeout-s", type=int, default=72_000)
    args = parser.parse_args()
    if args.task_id_start < 0 or args.target <= 0 or args.sleep_s <= 0:
        parser.error("task start must be nonnegative; target and sleep must be positive")
    if args.task_id_start + args.target > 100_000:
        parser.error("requested task interval exceeds the 100k manifest")
    if args.max_jobs_per_submit <= 0 or args.retry_timeout_s <= 0:
        parser.error("retry batch size and timeout must be positive")
    if not re.fullmatch(r"[0-9a-f]{64}", args.expected_tables_sha256):
        parser.error("--expected-tables-sha256 must be a lowercase SHA256")
    for name in ("expected_alpha", "expected_broadening_k"):
        value = getattr(args, name)
        if not math.isfinite(value) or value < 0:
            parser.error(f"--{name.replace('_', '-')} must be finite and nonnegative")
    return args


def load_assignments(path: Path) -> dict[int, dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 100_000:
        raise ValueError(f"{path}: expected 100000 tasks, found {len(rows)}")
    assignments = {int(row["task_id"]): row for row in rows}
    if set(assignments) != set(range(100_000)):
        raise ValueError(f"{path}: task IDs must be exactly 0..99999")
    return assignments


def float_matches(raw: str | None, expected: float) -> bool:
    try:
        value = float(raw) if raw is not None else math.nan
    except ValueError:
        return False
    return math.isclose(value, expected, rel_tol=0.0, abs_tol=1e-12)


def row_matches(
    row: dict[str, str],
    assignment: dict[str, str],
    *,
    variant: str,
    use_prehydro: bool,
    alpha: float,
    broadening_k: float,
    tables_sha256: str,
) -> bool:
    expected = {
        "task_id": assignment["task_id"],
        "seed": assignment["hard_seed"],
        "centrality": "C0-5",
        "hydro_slot": assignment["hydro_slot"],
        "hydro_event_id": assignment["hydro_event_id"],
        "hydro_ncoll": assignment["hydro_ncoll"],
        "hydro_payload_sha256": assignment["hydro_payload_sha256"],
        "variant": variant,
        "use_prehydro": str(int(use_prehydro)),
        "do_moliere": "1",
        "hadro_type": "1",
        "moliere_tables_sha256": tables_sha256,
        "returncode": "0",
        "timeout": "0",
    }
    return (
        all(row.get(key) == value for key, value in expected.items())
        and float_matches(row.get("energy_loss_alpha"), alpha)
        and float_matches(row.get("broadening_k"), broadening_k)
    )


def status_matches(
    status: dict[str, str],
    assignment: dict[str, str],
    *,
    alpha: float,
    broadening_k: float,
    tables_sha256: str,
) -> bool:
    expected = {
        "kind": "aa",
        "task_id": assignment["task_id"],
        "status": "success",
        "exit_code": "0",
        "run_prehydro_pair": "true",
        "run_prehydro_only": "false",
        "do_moliere": "true",
        "moliere_mode": "legacy_resolved",
        "moliere_tables_sha256": tables_sha256,
        "hydro_slot": assignment["hydro_slot"],
        "hydro_event_id": assignment["hydro_event_id"],
        "hydro_ncoll": assignment["hydro_ncoll"],
        "hydro_dir": assignment["hydro_dir"],
        "hydro_payload_key": assignment["hydro_payload_key"],
        "hydro_payload_sha256": assignment["hydro_payload_sha256"],
    }
    return (
        all(status.get(key) == value for key, value in expected.items())
        and float_matches(status.get("no_prehydro_alpha"), alpha)
        and float_matches(status.get("prehydro_alpha"), alpha)
        and float_matches(status.get("broadening_k"), broadening_k)
    )


def parse_input(archive: tarfile.TarFile, member_name: str) -> dict[str, str]:
    member = archive.getmember(member_name)
    stream = archive.extractfile(member)
    if stream is None:
        raise ValueError(f"archive member is not a regular file: {member_name}")
    values: dict[str, str] = {}
    for payload in stream.read().decode("utf-8", errors="strict").splitlines():
        if "=" in payload:
            key, value = payload.split("=", 1)
            values[key.strip()] = value.strip()
    return values


def input_matches(values: dict[str, str], *, use_prehydro: bool, alpha: float) -> bool:
    expected = {
        "do_quench": "true",
        "do_wake": "true",
        "do_elastic": "true",
        "do_lres": "false",
        "do_Moliere_on_unresolved_partons": "false",
        "do_Moliere_dynamic_unresolved_resolution": "false",
        "do_Moliere_dynamic_daughter_unresolved_resolution": "false",
        "do_Moliere_recursive_unresolved_resolution": "false",
        "hadro_type": "1",
        "compat_moliere_legacy_hydro": "true",
        "use_prehydro": "true" if use_prehydro else "false",
    }
    return (
        all(values.get(key) == value for key, value in expected.items())
        and bool(values.get("tables_path"))
        and float_matches(values.get("alpha"), alpha)
    )


def archive_is_complete(
    path: Path,
    assignment: dict[str, str],
    *,
    alpha: float,
    broadening_k: float,
    tables_sha256: str,
) -> bool:
    task_id = int(assignment["task_id"])
    task = f"task_{task_id:05d}"
    required_suffixes = {
        "baseline_summary": f"/{task}/summary.tsv",
        "baseline_hadrons": f"/{task}/HYBRID_Hadrons.out",
        "baseline_input": f"/{task}/hybrid_input.dat",
        "prehydro_summary": f"/{task}_prehydro/summary.tsv",
        "prehydro_hadrons": f"/{task}_prehydro/HYBRID_Hadrons.out",
        "prehydro_input": f"/{task}_prehydro/hybrid_input.dat",
        "pair_summary": f"/{task}_pair_summary.tsv",
    }
    try:
        with tarfile.open(path, "r:gz") as archive:
            matched: dict[str, tarfile.TarInfo] = {}
            members = archive.getmembers()
            for key, suffix in required_suffixes.items():
                candidates = [member for member in members if member.name.endswith(suffix)]
                if len(candidates) != 1 or not candidates[0].isfile():
                    return False
                matched[key] = candidates[0]
            if matched["baseline_hadrons"].size <= 0 or matched["prehydro_hadrons"].size <= 0:
                return False

            pair_rows = parse_tsv_member(archive, matched["pair_summary"].name)
            if len(pair_rows) != 2:
                return False
            rows_by_variant = {row.get("variant"): row for row in pair_rows}
            for variant, use_prehydro in (("no_prehydro", False), ("with_prehydro", True)):
                row = rows_by_variant.get(variant)
                if row is None or not row_matches(
                    row,
                    assignment,
                    variant=variant,
                    use_prehydro=use_prehydro,
                    alpha=alpha,
                    broadening_k=broadening_k,
                    tables_sha256=tables_sha256,
                ):
                    return False

            for key, variant, use_prehydro in (
                ("baseline_summary", "no_prehydro", False),
                ("prehydro_summary", "with_prehydro", True),
            ):
                rows = parse_tsv_member(archive, matched[key].name)
                if len(rows) != 1 or not row_matches(
                    rows[0],
                    assignment,
                    variant=variant,
                    use_prehydro=use_prehydro,
                    alpha=alpha,
                    broadening_k=broadening_k,
                    tables_sha256=tables_sha256,
                ):
                    return False

            if not input_matches(
                parse_input(archive, matched["baseline_input"].name),
                use_prehydro=False,
                alpha=alpha,
            ) or not input_matches(
                parse_input(archive, matched["prehydro_input"].name),
                use_prehydro=True,
                alpha=alpha,
            ):
                return False
    except (KeyError, OSError, tarfile.TarError, UnicodeError, ValueError):
        return False
    return True


def ready_ids(
    local_eos: Path,
    assignments: dict[int, dict[str, str]],
    expected_ids: set[int],
    *,
    alpha: float,
    broadening_k: float,
    tables_sha256: str,
) -> set[int]:
    ready: set[int] = set()
    for task_id in expected_ids:
        assignment = assignments[task_id]
        status_path = local_eos / f"status/aa/chunk_{task_id}.txt"
        output_path = local_eos / f"outputs/aa/chunk_{task_id}.tar.gz"
        if not status_path.is_file() or not output_path.is_file():
            continue
        if not status_matches(
            parse_status(status_path),
            assignment,
            alpha=alpha,
            broadening_k=broadening_k,
            tables_sha256=tables_sha256,
        ):
            continue
        if archive_is_complete(
            output_path,
            assignment,
            alpha=alpha,
            broadening_k=broadening_k,
            tables_sha256=tables_sha256,
        ):
            ready.add(task_id)
    return ready


def enforce_job_priority(args: argparse.Namespace) -> None:
    if args.dry_run:
        return
    constraint = f"{campaign_constraint(args.campaign)} && JobPrio != {args.job_priority}"
    query_command = (
        f"condor_q -name {args.schedd} -constraint '{constraint}' -af ClusterId"
    )
    query = subprocess.run(
        [args.cernctl, "run", "bash", "-lc", query_command],
        check=True,
        text=True,
        capture_output=True,
    )
    if not any(line.strip() for line in query.stdout.splitlines()):
        return
    command = (
        f"condor_qedit -name {args.schedd} -constraint '{constraint}' "
        f"JobPrio {args.job_priority}"
    )
    subprocess.run(
        [args.cernctl, "run", "bash", "-lc", command],
        check=True,
        text=True,
        capture_output=True,
    )


def supervise_once(args: argparse.Namespace) -> bool:
    assignments = load_assignments(args.task_manifest)
    expected_ids = set(range(args.task_id_start, args.task_id_start + args.target))
    local_eos = args.work / "eos_snapshot"
    seen_path = args.work / "moliere_pair_queue_seen.txt"
    states, active_aa, held_aa = query_jobs(args)
    enforce_job_priority(args)
    if active_aa and not seen_path.exists():
        seen_path.write_text(
            f"date={datetime.now().astimezone().isoformat(timespec='seconds')}\n"
        )
    sync_eos(args, local_eos)
    ready = ready_ids(
        local_eos,
        assignments,
        expected_ids,
        alpha=args.expected_alpha,
        broadening_k=args.expected_broadening_k,
        tables_sha256=args.expected_tables_sha256,
    )
    log(
        f"Condor states={dict(sorted(states.items()))}; "
        f"active_AA={len(active_aa)} strict_ready={len(ready)}/{args.target}"
    )
    active_aa -= remove_terminal_holds(args, held_aa, ready)
    if active_aa or not seen_path.exists():
        return False
    missing = expected_ids - ready
    if missing:
        submit_retry(args, missing)
        return False

    marker = args.work / "moliere_pair_strict_complete.txt"
    marker.write_text(
        f"date={datetime.now().astimezone().isoformat(timespec='seconds')}\n"
        f"accepted_pairs={len(ready)}\n"
        f"alpha={args.expected_alpha}\n"
        f"moliere_tables_sha256={args.expected_tables_sha256}\n"
    )
    log("paired Moliere strict target complete")
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
