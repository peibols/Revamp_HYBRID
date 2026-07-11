#!/usr/bin/env python3
"""Drain, audit, and repair the guarded OO paired-AA campaigns to 50k."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import re
import subprocess
import time


DEFAULT_ROOT = Path("/raid5/data/yjlee/hybrid_dev")
RECOVERY_CAMPAIGN = (
    "hybrid_oo5360_c0_5_no_moliere_paired_public2509_planB_"
    "17369AA_recovery_guard_20260710"
)
EXTENSION_CAMPAIGN = (
    "hybrid_oo5360_c0_5_no_moliere_paired_public2509_planB_"
    "20kAA_extension_to50k_guard_20260711"
)
ORIGINAL_WORK = "test/oo5360_planb_public_10k_20260710"
RETAINED_WORK = "test/oo5360_planb_public_20k_continuation_to30k_20260710"
RECOVERY_WORK = "test/oo5360_planb_public_20k_recovery_guard_20260710"
EXTENSION_WORK = "test/oo5360_planb_public_20k_extension_to50k_20260711"
PP_WORK = "test/tmp_oo_10k_prehydro_raa_20260709/local_eos"


@dataclass(frozen=True)
class LiveRoot:
    label: str
    campaign: str
    work: Path
    local_eos: Path
    expected_id_files: tuple[Path, ...]
    submit_template: Path

    @property
    def eos_base(self) -> str:
        return f"/eos/user/y/yjlee/{self.campaign}"

    @property
    def afs_work(self) -> str:
        return f"/afs/cern.ch/user/y/yjlee/cernLxplus_jobs/{self.campaign}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--target", type=int, default=50000)
    parser.add_argument("--sleep-s", type=int, default=900)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--cernctl", default="/data/yjlee/cernLxplus/cernctl")
    parser.add_argument("--cern-remote", default="lxplus")
    parser.add_argument("--schedd", default="bigbird103.cern.ch")
    args = parser.parse_args()
    if args.target <= 0:
        parser.error("--target must be positive")
    if args.sleep_s <= 0:
        parser.error("--sleep-s must be positive")
    return args


def read_ids(paths: tuple[Path, ...]) -> set[int]:
    ids: set[int] = set()
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        for raw in path.read_text().splitlines():
            value = raw.strip()
            if value:
                ids.add(int(value))
    return ids


def successful_ids(local_eos: Path) -> set[int]:
    status_dir = local_eos / "status/aa"
    ids: set[int] = set()
    if not status_dir.exists():
        return ids
    for path in status_dir.glob("chunk_*.txt"):
        if "status=success" not in path.read_text(errors="replace").splitlines():
            continue
        match = re.fullmatch(r"chunk_(\d+)\.txt", path.name)
        if match:
            ids.add(int(match.group(1)))
    return ids


def output_ids(local_eos: Path) -> set[int]:
    output_dir = local_eos / "outputs/aa"
    ids: set[int] = set()
    if not output_dir.exists():
        return ids
    for path in output_dir.glob("chunk_*.tar.gz"):
        match = re.fullmatch(r"chunk_(\d+)\.tar\.gz", path.name)
        if match and path.stat().st_size > 0:
            ids.add(int(match.group(1)))
    return ids


def ready_ids(local_eos: Path) -> set[int]:
    return successful_ids(local_eos) & output_ids(local_eos)


def timestamp() -> str:
    return datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z")


def log(message: str) -> None:
    print(f"{datetime.now().astimezone().isoformat(timespec='seconds')} {message}", flush=True)


def active_job_states(
    root: LiveRoot, *, cernctl: str, schedd: str
) -> Counter[int]:
    constraint = f'regexp("{root.campaign}", Cmd)'
    result = subprocess.run(
        [
            cernctl,
            "run",
            "bash",
            "-lc",
            f"condor_q -name {schedd} -constraint '{constraint}' -af JobStatus",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    return Counter(int(line) for line in result.stdout.splitlines() if line.strip())


def sync_root(root: LiveRoot, cern_remote: str) -> None:
    for member in ("status", "outputs"):
        destination = root.local_eos / member
        destination.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                "rsync",
                "-a",
                "--partial",
                f"{cern_remote}:{root.eos_base}/{member}/",
                f"{destination}/",
            ],
            check=True,
        )


def render_retry_submit(template: str, id_file_name: str, tag: str) -> str:
    queue_pattern = re.compile(r"^queue chunk_id from .+$", re.MULTILINE)
    if not queue_pattern.search(template):
        raise ValueError("submit template has no chunk-id queue statement")
    rendered = queue_pattern.sub(f"queue chunk_id from {id_file_name}", template)
    for field in ("output", "error", "log"):
        pattern = re.compile(rf"^({field}\s*=\s*log/)([^\n]+)$", re.MULTILINE)
        match = pattern.search(rendered)
        if not match:
            raise ValueError(f"submit template has no {field} log path")
        old = match.group(2)
        rendered = pattern.sub(rf"\g<1>reconcile_{tag}_{old}", rendered, count=1)
    return rendered


def append_submission(
    ledger: Path,
    *,
    stamp: str,
    root: LiveRoot,
    ids: set[int],
    cluster: str,
) -> None:
    ledger.parent.mkdir(parents=True, exist_ok=True)
    exists = ledger.exists()
    with ledger.open("a", newline="") as handle:
        fieldnames = [
            "date",
            "root",
            "cluster",
            "job_count",
            "task_min",
            "task_max",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        if not exists:
            writer.writeheader()
        writer.writerow(
            {
                "date": stamp,
                "root": root.label,
                "cluster": cluster,
                "job_count": len(ids),
                "task_min": min(ids),
                "task_max": max(ids),
            }
        )


def submit_ids(
    root: LiveRoot,
    ids: set[int],
    *,
    extension_work: Path,
    cernctl: str,
    cern_remote: str,
    dry_run: bool,
) -> str:
    if not ids:
        raise ValueError("cannot submit an empty retry set")
    stamp = timestamp()
    reconcile_dir = extension_work / "reconcile"
    reconcile_dir.mkdir(parents=True, exist_ok=True)
    id_path = reconcile_dir / f"{root.label}_{stamp}_ids.txt"
    submit_path = reconcile_dir / f"{root.label}_{stamp}.sub"
    id_path.write_text("".join(f"{task_id}\n" for task_id in sorted(ids)))
    submit_path.write_text(
        render_retry_submit(root.submit_template.read_text(), id_path.name, stamp)
    )
    if dry_run:
        log(f"dry run: would submit {len(ids)} {root.label} task(s)")
        return "dry-run"
    subprocess.run(
        [
            "scp",
            "-q",
            "-o",
            "BatchMode=yes",
            str(id_path),
            str(submit_path),
            f"{cern_remote}:{root.afs_work}/",
        ],
        check=True,
    )
    command = (
        "source /etc/profile.d/modules.sh 2>/dev/null || true; "
        "module load lxbatch/eossubmit >/dev/null 2>&1; "
        f"cd {root.afs_work} && condor_submit {submit_path.name}"
    )
    result = subprocess.run(
        [cernctl, "run", "bash", "-lc", command],
        check=True,
        text=True,
        capture_output=True,
    )
    match = re.search(r"submitted to cluster (\d+)", result.stdout)
    if not match:
        raise RuntimeError(f"could not parse retry cluster from: {result.stdout}")
    cluster = match.group(1)
    append_submission(
        extension_work / "reconcile_submissions.tsv",
        stamp=stamp,
        root=root,
        ids=ids,
        cluster=cluster,
    )
    log(f"submitted {len(ids)} {root.label} task(s) as cluster {cluster}")
    return cluster


def run_strict_analysis(
    *,
    root: Path,
    extension_work: Path,
    recovery: LiveRoot,
    extension: LiveRoot,
) -> tuple[int, dict[Path, set[int]], Path]:
    out_dir = extension_work / "reconcile_analysis" / timestamp()
    analyzer = extension_work / "analysis/analyze_oo_prehydro_pair.py"
    original = root / ORIGINAL_WORK / "eos_snapshot"
    retained = root / RETAINED_WORK / "eos_snapshot"
    pp = root / PP_WORK
    subprocess.run(
        [
            "python3",
            str(analyzer),
            "--local-eos",
            str(recovery.local_eos),
            "--additional-aa-local-eos",
            str(original),
            "--additional-aa-local-eos",
            str(retained),
            "--additional-aa-local-eos",
            str(extension.local_eos),
            "--pp-local-eos",
            str(pp),
            "--out-dir",
            str(out_dir),
            "--require-paired-aa",
            "--aa-events-per-chunk",
            "1",
        ],
        check=True,
    )
    with (out_dir / "aa_sources.tsv").open() as handle:
        accepted = sum(
            int(row["chunks_used"])
            for row in csv.DictReader(handle, delimiter="\t")
        )
    rejected: dict[Path, set[int]] = {}
    with (out_dir / "aa_rejections.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            source = Path(row["aa_local_eos"]).resolve()
            rejected.setdefault(source, set()).add(int(row["chunk_id"]))
    return accepted, rejected, out_dir


def build_roots(root: Path) -> tuple[LiveRoot, LiveRoot, Path]:
    recovery_work = root / RECOVERY_WORK
    extension_work = root / EXTENSION_WORK
    recovery = LiveRoot(
        label="recovery",
        campaign=RECOVERY_CAMPAIGN,
        work=recovery_work,
        local_eos=recovery_work / "eos_snapshot",
        expected_id_files=(
            recovery_work / "aa_chunk_ids.txt",
            recovery_work / "aa_chunk_id_3478.txt",
        ),
        submit_template=recovery_work / "oo_no_moliere_aa.sub",
    )
    extension = LiveRoot(
        label="extension",
        campaign=EXTENSION_CAMPAIGN,
        work=extension_work,
        local_eos=extension_work / "eos_snapshot",
        expected_id_files=(
            extension_work / "aa_chunk_ids.txt",
            extension_work / "aa_chunk_id_20000.txt",
            extension_work / "aa_chunk_ids_static_replacements_1788.txt",
        ),
        submit_template=extension_work / "oo_no_moliere_aa.sub",
    )
    return recovery, extension, extension_work


def supervise_once(args: argparse.Namespace) -> bool:
    recovery, extension, extension_work = build_roots(args.root)
    live_roots = (recovery, extension)
    active = False
    for live_root in live_roots:
        states = active_job_states(
            live_root, cernctl=args.cernctl, schedd=args.schedd
        )
        log(f"{live_root.label} Condor states: {dict(sorted(states.items()))}")
        active = active or bool(states)
    if active:
        return False

    missing_by_root: dict[LiveRoot, set[int]] = {}
    for live_root in live_roots:
        sync_root(live_root, args.cern_remote)
        expected = read_ids(live_root.expected_id_files)
        missing = expected - ready_ids(live_root.local_eos)
        log(
            f"{live_root.label} file audit: "
            f"ready={len(expected) - len(missing)}/{len(expected)}"
        )
        if missing:
            missing_by_root[live_root] = missing
    if missing_by_root:
        for live_root, ids in missing_by_root.items():
            submit_ids(
                live_root,
                ids,
                extension_work=extension_work,
                cernctl=args.cernctl,
                cern_remote=args.cern_remote,
                dry_run=args.dry_run,
            )
        return False

    accepted, rejected, out_dir = run_strict_analysis(
        root=args.root,
        extension_work=extension_work,
        recovery=recovery,
        extension=extension,
    )
    log(f"strict analysis accepted {accepted}/{args.target}; output={out_dir}")
    if accepted >= args.target:
        marker = extension_work / "combined_50k_strict_complete.txt"
        marker.write_text(
            f"date={datetime.now().astimezone().isoformat(timespec='seconds')}\n"
            f"accepted_pairs={accepted}\n"
            f"target_pairs={args.target}\n"
            f"analysis_dir={out_dir}\n"
        )
        return True

    live_rejected: dict[LiveRoot, set[int]] = {}
    for live_root in live_roots:
        ids = rejected.get(live_root.local_eos.resolve(), set())
        if ids:
            live_rejected[live_root] = ids
    if not live_rejected:
        raise RuntimeError(
            f"strict count is short by {args.target - accepted}, but no live-root "
            "rejection can be resubmitted"
        )
    for live_root, ids in live_rejected.items():
        submit_ids(
            live_root,
            ids,
            extension_work=extension_work,
            cernctl=args.cernctl,
            cern_remote=args.cern_remote,
            dry_run=args.dry_run,
        )
    return False


def main() -> int:
    args = parse_args()
    while True:
        try:
            complete = supervise_once(args)
        except (OSError, RuntimeError, subprocess.SubprocessError, ValueError) as error:
            log(f"supervisor error: {type(error).__name__}: {error}")
            complete = False
        if complete or args.once:
            return 0 if complete else 1
        time.sleep(args.sleep_s)


if __name__ == "__main__":
    raise SystemExit(main())
