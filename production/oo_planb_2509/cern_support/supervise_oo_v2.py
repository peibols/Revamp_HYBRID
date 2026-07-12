#!/usr/bin/env python3
"""Monitor, prefix-audit, and repair the 500-hydro OO v2 campaign."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime
import io
from pathlib import Path
import re
import subprocess
import tarfile
import time


DEFAULT_ROOT = Path("/raid5/data/yjlee/hybrid_dev")
DEFAULT_WORK = DEFAULT_ROOT / "test/oo5360_v2_500hydro_50k_20260711"
DEFAULT_SOURCE = DEFAULT_ROOT / "wt_main_moliere_lres_integration_clean"
DEFAULT_CAMPAIGN = (
    "hybrid_oo5360_c0_5_500hydro_no_moliere_paired_public2509_"
    "planB_v2_50kAA_20260711"
)
DEFAULT_AFS_WORK = Path(
    "/afs/cern.ch/user/y/yjlee/oo_v2_500hydro_50k_20260711"
)
DEFAULT_PP = DEFAULT_ROOT / "test/tmp_oo_10k_prehydro_raa_20260709/local_eos"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--campaign", default=DEFAULT_CAMPAIGN)
    parser.add_argument("--task-id-start", type=int, default=0)
    parser.add_argument("--target", type=int, default=50_000)
    parser.add_argument("--block-size", type=int, default=5_000)
    parser.add_argument("--sleep-s", type=int, default=900)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--cernctl", default="/data/yjlee/cernLxplus/cernctl")
    parser.add_argument("--cern-remote", default="lxplus")
    parser.add_argument("--schedd", default="bigbird101.cern.ch")
    parser.add_argument("--afs-work", type=Path, default=DEFAULT_AFS_WORK)
    parser.add_argument("--max-jobs-per-submit", type=int, default=10_000)
    parser.add_argument("--retry-timeout-s", type=int, default=72_000)
    parser.add_argument("--pp-local-eos", type=Path, default=DEFAULT_PP)
    args = parser.parse_args()
    if args.task_id_start < 0 or args.target <= 0 or args.block_size <= 0:
        parser.error("--task-id-start must be nonnegative; target/block size positive")
    if args.target % args.block_size:
        parser.error("--target must be divisible by --block-size")
    if (
        args.sleep_s <= 0
        or args.max_jobs_per_submit <= 0
        or args.retry_timeout_s <= 0
    ):
        parser.error(
            "--sleep-s, --max-jobs-per-submit, and --retry-timeout-s must be positive"
        )
    return args


def log(message: str) -> None:
    print(
        f"{datetime.now().astimezone().isoformat(timespec='seconds')} {message}",
        flush=True,
    )


def parse_status(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw in path.read_text(errors="replace").splitlines():
        if "=" in raw:
            key, value = raw.split("=", 1)
            values[key] = value
    return values


def successful_ids(local_eos: Path) -> set[int]:
    ids: set[int] = set()
    for path in (local_eos / "status/aa").glob("chunk_*.txt"):
        match = re.fullmatch(r"chunk_(\d+)\.txt", path.name)
        if match and parse_status(path).get("status") == "success":
            ids.add(int(match.group(1)))
    return ids


def parse_tsv_member(archive: tarfile.TarFile, name: str) -> list[dict[str, str]]:
    member = archive.getmember(name)
    stream = archive.extractfile(member)
    if stream is None:
        raise ValueError(f"archive member is not a regular file: {name}")
    with io.TextIOWrapper(stream, encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def paired_archive_is_complete(path: Path, task_id: int) -> bool:
    task = f"task_{task_id:05d}"
    required_suffixes = {
        "baseline_summary": f"/{task}/summary.tsv",
        "baseline_hadrons": f"/{task}/HYBRID_Hadrons.out",
        "prehydro_summary": f"/{task}_prehydro/summary.tsv",
        "prehydro_hadrons": f"/{task}_prehydro/HYBRID_Hadrons.out",
        "pair_summary": f"/{task}_pair_summary.tsv",
    }
    try:
        with tarfile.open(path, "r:gz") as archive:
            members = archive.getmembers()
            matched: dict[str, tarfile.TarInfo] = {}
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
            if {row.get("variant") for row in pair_rows} != {
                "no_prehydro",
                "with_prehydro",
            }:
                return False
            for row in pair_rows:
                if (
                    row.get("task_id") != str(task_id)
                    or row.get("returncode") != "0"
                    or row.get("timeout") != "0"
                ):
                    return False

            for key, expected_variant in (
                ("baseline_summary", "no_prehydro"),
                ("prehydro_summary", "with_prehydro"),
            ):
                rows = parse_tsv_member(archive, matched[key].name)
                if len(rows) != 1:
                    return False
                row = rows[0]
                if (
                    row.get("variant") != expected_variant
                    or row.get("task_id") != str(task_id)
                    or row.get("returncode") != "0"
                    or row.get("timeout") != "0"
                ):
                    return False
    except (KeyError, OSError, tarfile.TarError, UnicodeError, ValueError):
        return False
    return True


def output_ids(local_eos: Path) -> set[int]:
    ids: set[int] = set()
    for path in (local_eos / "outputs/aa").glob("chunk_*.tar.gz"):
        match = re.fullmatch(r"chunk_(\d+)\.tar\.gz", path.name)
        if not match or path.stat().st_size <= 0:
            continue
        task_id = int(match.group(1))
        if paired_archive_is_complete(path, task_id):
            ids.add(task_id)
    return ids


def ready_ids(local_eos: Path) -> set[int]:
    return successful_ids(local_eos) & output_ids(local_eos)


def query_jobs(
    args: argparse.Namespace,
) -> tuple[Counter[int], set[int], set[int]]:
    constraint = f'regexp("{args.campaign}", Environment)'
    command = (
        f"condor_q -name {args.schedd} -constraint '{constraint}' "
        "-af JobStatus Args"
    )
    result = subprocess.run(
        [args.cernctl, "run", "bash", "-lc", command],
        check=True,
        text=True,
        capture_output=True,
    )
    states: Counter[int] = Counter()
    active_aa: set[int] = set()
    held_aa: set[int] = set()
    for raw in result.stdout.splitlines():
        fields = raw.split(maxsplit=1)
        if len(fields) != 2:
            continue
        status = int(fields[0])
        states[status] += 1
        if not fields[1].isdigit():
            continue
        task_id = int(fields[1])
        if status in {1, 2, 5}:
            active_aa.add(task_id)
        if status == 5:
            held_aa.add(task_id)
    return states, active_aa, held_aa


def remove_terminal_holds(
    args: argparse.Namespace,
    held_aa: set[int],
    ready: set[int],
) -> set[int]:
    if not held_aa:
        return set()
    strict_ready = held_aa & ready
    needs_retry = held_aa - ready
    if args.dry_run:
        log(
            f"dry run: would remove {len(held_aa)} terminal held task(s): "
            f"{len(strict_ready)} strict-ready, {len(needs_retry)} need retry"
        )
        return set()
    constraint = (
        f'regexp("{args.campaign}", Environment) && JobStatus == 5'
    )
    command = (
        f"condor_rm -name {args.schedd} -constraint '{constraint}'"
    )
    subprocess.run(
        [args.cernctl, "run", "bash", "-lc", command],
        check=True,
        text=True,
        capture_output=True,
    )
    log(
        f"removed {len(held_aa)} terminal held task(s): "
        f"{len(strict_ready)} strict-ready, {len(needs_retry)} need retry"
    )
    return set(held_aa)


def sync_eos(args: argparse.Namespace, local_eos: Path) -> None:
    eos_base = f"/eos/user/y/yjlee/{args.campaign}"
    for member in ("status", "outputs"):
        destination = local_eos / member
        destination.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                "rsync",
                "-a",
                "--partial",
                f"{args.cern_remote}:{eos_base}/{member}/",
                f"{destination}/",
            ],
            check=True,
        )


def completed_milestones(path: Path) -> set[int]:
    if not path.is_file():
        return set()
    with path.open() as handle:
        return {
            int(row["task_limit"])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def append_milestone(
    path: Path,
    task_limit: int,
    analysis_dir: Path,
    *,
    task_start: int,
    target: int,
) -> None:
    exists = path.exists()
    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("date", "task_limit", "percent", "analysis_dir"),
            delimiter="\t",
            lineterminator="\n",
        )
        if not exists:
            writer.writeheader()
        writer.writerow(
            {
                "date": datetime.now().astimezone().isoformat(timespec="seconds"),
                "task_limit": task_limit,
                "percent": (task_limit - task_start) * 100 // target,
                "analysis_dir": analysis_dir,
            }
        )


def run_analysis(args: argparse.Namespace, local_eos: Path, task_limit: int) -> Path:
    out_dir = (
        args.work
        / "analysis_v2_milestones"
        / f"tasks_{args.task_id_start:05d}_{task_limit:05d}"
    )
    analyzer = (
        args.source
        / "production/oo_planb_2509/analysis/analyze_oo_prehydro_pair.py"
    )
    manifest = args.work / "hydro_prepared/aa_task_manifest.tsv"
    subprocess.run(
        [
            "python3",
            str(analyzer),
            "--local-eos",
            str(local_eos),
            "--pp-local-eos",
            str(args.pp_local_eos),
            "--out-dir",
            str(out_dir),
            "--require-paired-aa",
            "--aa-events-per-chunk",
            "1",
            "--aa-task-manifest",
            str(manifest),
            "--aa-task-start",
            str(args.task_id_start),
            "--aa-task-limit",
            str(task_limit),
            "--require-complete-aa-prefix",
        ],
        check=True,
    )
    return out_dir


def render_retry_submit(
    template: str,
    id_file: str,
    tag: str,
    retry_timeout_s: int | None = None,
) -> str:
    queue = re.compile(r"^queue chunk_id from .+$", re.MULTILINE)
    if not queue.search(template):
        raise ValueError("AA submit template has no chunk-id queue statement")
    rendered = queue.sub(f"queue chunk_id from {id_file}", template)
    transfer_pattern = re.compile(
        r"^transfer_output_files\s*=.*$", re.MULTILINE
    )
    if transfer_pattern.search(rendered):
        rendered = transfer_pattern.sub(
            'transfer_output_files = ""', rendered, count=1
        )
    else:
        transfer_marker = re.compile(
            r"^(when_to_transfer_output\s*=.*)$", re.MULTILINE
        )
        if not transfer_marker.search(rendered):
            raise ValueError(
                "AA submit template has no output-transfer marker"
            )
        rendered = transfer_marker.sub(
            '\\1\ntransfer_output_files = ""', rendered, count=1
        )
    for field in ("output", "error"):
        pattern = re.compile(rf"^({field}\s*=\s*log/)([^\n]+)$", re.MULTILINE)
        if pattern.search(rendered):
            rendered = pattern.sub(
                rf"\g<1>v2_retry_{tag}_\2", rendered, count=1
            )
            continue
        null_pattern = re.compile(
            rf"^{field}\s*=\s*/dev/null\s*$", re.MULTILINE
        )
        if not null_pattern.search(rendered):
            raise ValueError(
                f"AA submit template has no supported {field} path"
            )
    log_pattern = re.compile(
        r"^log\s*=\s*(?:log/[^\n]+|/dev/null)\s*$", re.MULTILINE
    )
    if not log_pattern.search(rendered):
        raise ValueError("AA submit template has no supported log path")
    rendered = log_pattern.sub("log = /dev/null", rendered, count=1)
    if retry_timeout_s is not None:
        timeout_pattern = re.compile(r"\bTIMEOUT_S=\d+\b")
        if len(timeout_pattern.findall(rendered)) != 1:
            raise ValueError("AA submit template must define TIMEOUT_S exactly once")
        rendered = timeout_pattern.sub(
            f"TIMEOUT_S={retry_timeout_s}", rendered, count=1
        )
    return rendered


def batched_ids(ids: set[int], batch_size: int) -> list[list[int]]:
    ordered = sorted(ids)
    return [
        ordered[start : start + batch_size]
        for start in range(0, len(ordered), batch_size)
    ]


def submit_retry(args: argparse.Namespace, ids: set[int]) -> None:
    if not ids:
        return
    stamp = datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z")
    retry_dir = args.work / "v2_retries"
    retry_dir.mkdir(parents=True, exist_ok=True)
    template = (args.work / "oo_no_moliere_aa.sub").read_text()
    batches = batched_ids(ids, args.max_jobs_per_submit)
    if args.dry_run:
        log(
            f"dry run: would resubmit {len(ids)} task(s) "
            f"in {len(batches)} cluster(s)"
        )
        return
    for part, batch in enumerate(batches):
        tag = f"{stamp}_part{part:02d}"
        id_path = retry_dir / f"retry_{tag}_ids.txt"
        submit_path = retry_dir / f"retry_{tag}.sub"
        id_path.write_text("".join(f"{task_id}\n" for task_id in batch))
        submit_path.write_text(
            render_retry_submit(
                template,
                id_path.name,
                tag,
                retry_timeout_s=args.retry_timeout_s,
            )
        )
        subprocess.run(
            [
                "scp",
                "-q",
                "-o",
                "BatchMode=yes",
                str(id_path),
                str(submit_path),
                f"{args.cern_remote}:{args.afs_work}/",
            ],
            check=True,
        )
        command = (
            "source /etc/profile.d/modules.sh 2>/dev/null || true; "
            "module load lxbatch/eossubmit >/dev/null 2>&1; "
            f"cd {args.afs_work} && "
            f"condor_submit -name {args.schedd} {submit_path.name}"
        )
        result = subprocess.run(
            [args.cernctl, "run", "bash", "-lc", command],
            check=True,
            text=True,
            capture_output=True,
        )
        match = re.search(r"submitted to cluster (\d+)", result.stdout)
        if not match:
            raise RuntimeError(f"could not parse retry cluster: {result.stdout}")
        with (retry_dir / "retry_submissions.tsv").open("a") as handle:
            handle.write(
                f"{datetime.now().astimezone().isoformat(timespec='seconds')}\t"
                f"{match.group(1)}\t{len(batch)}\t{min(batch)}\t{max(batch)}\t"
                f"{args.retry_timeout_s}\n"
            )
        log(
            f"resubmitted {len(batch)} task(s) as cluster {match.group(1)} "
            f"({part + 1}/{len(batches)}), timeout={args.retry_timeout_s}s"
        )


def rejected_ids(analysis_dir: Path) -> set[int]:
    path = analysis_dir / "aa_rejections.tsv"
    if not path.is_file():
        return set()
    with path.open() as handle:
        return {
            int(row["chunk_id"])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def supervise_once(args: argparse.Namespace) -> bool:
    local_eos = args.work / "eos_snapshot"
    seen_path = args.work / "v2_aa_queue_seen.txt"
    state_path = args.work / "v2_milestones.tsv"
    states, active_aa, held_aa = query_jobs(args)
    if active_aa and not seen_path.exists():
        seen_path.write_text(
            f"date={datetime.now().astimezone().isoformat(timespec='seconds')}\n"
        )
    sync_eos(args, local_eos)
    task_stop = args.task_id_start + args.target
    expected_ids = set(range(args.task_id_start, task_stop))
    ready = ready_ids(local_eos) & expected_ids
    log(
        f"Condor states={dict(sorted(states.items()))}; "
        f"active_AA={len(active_aa)} ready={len(ready)}/{args.target}"
    )
    active_aa -= remove_terminal_holds(args, held_aa, ready)

    done = completed_milestones(state_path)
    for completed_count in range(
        args.block_size, args.target + 1, args.block_size
    ):
        task_limit = args.task_id_start + completed_count
        if task_limit in done:
            continue
        if not set(range(args.task_id_start, task_limit)).issubset(ready):
            break
        try:
            out_dir = run_analysis(args, local_eos, task_limit)
        except subprocess.CalledProcessError as error:
            log(f"strict prefix analysis failed for {task_limit}: {error}")
            break
        append_milestone(
            state_path,
            task_limit,
            out_dir,
            task_start=args.task_id_start,
            target=args.target,
        )
        log(f"strict milestone complete: {completed_count}/{args.target}")

    if active_aa or not seen_path.exists():
        return False

    missing = expected_ids - ready
    if missing:
        submit_retry(args, missing)
        return False

    final_dir = (
        args.work
        / "analysis_v2_milestones"
        / f"tasks_{args.task_id_start:05d}_{task_stop:05d}"
    )
    try:
        if task_stop not in completed_milestones(state_path):
            final_dir = run_analysis(args, local_eos, task_stop)
            append_milestone(
                state_path,
                task_stop,
                final_dir,
                task_start=args.task_id_start,
                target=args.target,
            )
    except subprocess.CalledProcessError:
        malformed = rejected_ids(final_dir)
        if not malformed:
            raise RuntimeError("full strict audit failed without resubmittable task IDs")
        submit_retry(args, malformed)
        return False

    (args.work / "v2_50k_strict_complete.txt").write_text(
        f"date={datetime.now().astimezone().isoformat(timespec='seconds')}\n"
        f"task_id_start={args.task_id_start}\n"
        f"task_id_stop={task_stop}\n"
        f"accepted_pairs={args.target}\nanalysis_dir={final_dir}\n"
    )
    log("v2 strict target complete")
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
