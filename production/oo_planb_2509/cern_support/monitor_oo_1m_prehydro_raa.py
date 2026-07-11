#!/usr/bin/env python3
"""Monitor the OO 1M paired prehydro campaign and update the Overleaf deck."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
from pathlib import Path
import re
import shutil
import subprocess
import time


DEFAULT_ROOT = Path("/raid5/data/yjlee/hybrid_dev")
DEFAULT_WORK = DEFAULT_ROOT / "test/oo5360_no_moliere_raa_20260606"
DEFAULT_OVERLEAF = DEFAULT_ROOT / "overleaf_69d5314739aed083fc3bbb0a"
DEFAULT_CAMPAIGN = "hybrid_oo5360_c0_5_no_moliere_paired_prehydro_100kAA_1mPP_aa1evt_20260710"
DEFAULT_FIGURE = "report/figures/20260709-oo-c0-5-paired-prehydro-raa-current.pdf"
BEGIN_MARKER = "% OO_AUTO_STATUS_BEGIN"
END_MARKER = "% OO_AUTO_STATUS_END"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--overleaf", type=Path, default=DEFAULT_OVERLEAF)
    parser.add_argument("--campaign", default=DEFAULT_CAMPAIGN)
    parser.add_argument("--eos-base")
    parser.add_argument("--local-eos", type=Path)
    parser.add_argument(
        "--additional-aa-local-eos",
        action="append",
        type=Path,
        default=[],
        help="additional local AA snapshot included in combined milestones",
    )
    parser.add_argument(
        "--additional-aa-chunks",
        action="append",
        type=int,
        default=[],
        help="expected chunks for each --additional-aa-local-eos",
    )
    parser.add_argument(
        "--sync-additional-aa",
        action="append",
        nargs=2,
        metavar=("EOS_BASE", "LOCAL_EOS"),
        default=[],
        help=(
            "additional live AA EOS root and its matching local snapshot; "
            "requires --sync-via-cernctl"
        ),
    )
    parser.add_argument(
        "--pp-local-eos",
        type=Path,
        help="optional retained snapshot containing a separate pp denominator",
    )
    parser.add_argument("--analysis-out-root", type=Path)
    parser.add_argument("--tex", type=Path, default=Path("report/20260709-OO.tex"))
    parser.add_argument("--aa-chunks", type=int, default=100000)
    parser.add_argument("--aa-events", type=int, default=1)
    parser.add_argument("--pp-chunks", type=int, default=1000)
    parser.add_argument("--pp-events", type=int, default=1000)
    parser.add_argument("--milestones", default="10,20,30,40,50,60,70,80,90,100")
    parser.add_argument("--sleep-s", type=int, default=900)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--copy", action="store_true", help="copy the EOS tree to --local-eos before checking status")
    parser.add_argument("--sync-via-cernctl", action="store_true", help="sync EOS snapshots through lxplus using cernctl and scp")
    parser.add_argument("--cernctl", default="/data/yjlee/cernLxplus/cernctl")
    parser.add_argument("--cern-remote", default="lxplus")
    parser.add_argument("--remote-stage", default="/tmp/yjlee_oo_1m_monitor_snapshot.tar.gz")
    parser.add_argument("--renew-kerberos", action="store_true")
    parser.add_argument("--commit-overleaf", action="store_true")
    parser.add_argument("--push-overleaf", action="store_true")
    args = parser.parse_args()
    if args.aa_events != 1:
        parser.error("--aa-events must be 1 for paired OO production")
    if len(args.additional_aa_local_eos) != len(args.additional_aa_chunks):
        parser.error(
            "--additional-aa-local-eos and --additional-aa-chunks "
            "must be repeated the same number of times"
        )
    if any(value <= 0 for value in args.additional_aa_chunks):
        parser.error("--additional-aa-chunks values must be positive")
    args.sync_additional_aa = [
        (eos_base, Path(local_eos))
        for eos_base, local_eos in args.sync_additional_aa
    ]
    if args.sync_additional_aa and not args.sync_via_cernctl:
        parser.error("--sync-additional-aa requires --sync-via-cernctl")
    sync_paths = [local_eos for _, local_eos in args.sync_additional_aa]
    if len(sync_paths) != len(set(sync_paths)):
        parser.error("each --sync-additional-aa local snapshot must be unique")
    unknown_sync_paths = set(sync_paths) - set(args.additional_aa_local_eos)
    if unknown_sync_paths:
        parser.error(
            "each --sync-additional-aa LOCAL_EOS must also be supplied with "
            "--additional-aa-local-eos"
        )
    return args


def successful_chunks(local_eos: Path, kind: str) -> set[int]:
    status_dir = local_eos / "status" / kind
    chunks: set[int] = set()
    if not status_dir.exists():
        return chunks
    for path in status_dir.glob("chunk_*.txt"):
        text = path.read_text(errors="replace")
        if "status=success" not in text:
            continue
        match = re.search(r"chunk_(\d+)\.txt$", path.name)
        if match:
            chunks.add(int(match.group(1)))
    return chunks


def available_chunks(local_eos: Path, kind: str) -> set[int]:
    """Return successful chunks whose output archive is present locally."""
    output_dir = local_eos / "outputs" / kind
    outputs: set[int] = set()
    if output_dir.exists():
        for path in output_dir.glob("chunk_*.tar.gz"):
            match = re.fullmatch(r"chunk_(\d+)\.tar\.gz", path.name)
            if match:
                outputs.add(int(match.group(1)))
    return successful_chunks(local_eos, kind) & outputs


def copy_from_eos(eos_base: str, local_eos: Path) -> None:
    local_eos.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["xrdcp", "-r", "-f", f"root://eosuser.cern.ch/{eos_base.rstrip('/')}/", str(local_eos)],
        check=True,
    )


def sync_via_cernctl(
    *,
    cernctl: str,
    cern_remote: str,
    eos_base: str,
    local_eos: Path,
    remote_stage: str,
    include_outputs: bool,
    renew_kerberos: bool,
) -> None:
    local_eos.mkdir(parents=True, exist_ok=True)
    if renew_kerberos:
        subprocess.run([cernctl, "kerberos"], check=True)
    members = ["status"]
    if include_outputs:
        members.append("outputs")
    # Keep the legacy argument so existing handoff commands remain compatible.
    del remote_stage
    for member in members:
        destination = local_eos / member
        destination.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            [
                "rsync",
                "-a",
                "--partial",
                f"{cern_remote}:{eos_base.rstrip('/')}/{member}/",
                f"{destination}/",
            ],
            check=True,
        )


def sync_configured_aa_sources(
    args: argparse.Namespace,
    *,
    eos_base: str,
    local_eos: Path,
    include_outputs: bool,
) -> None:
    sync_via_cernctl(
        cernctl=args.cernctl,
        cern_remote=args.cern_remote,
        eos_base=eos_base,
        local_eos=local_eos,
        remote_stage=args.remote_stage,
        include_outputs=include_outputs,
        renew_kerberos=args.renew_kerberos,
    )
    for index, (additional_eos_base, additional_local_eos) in enumerate(
        args.sync_additional_aa
    ):
        sync_via_cernctl(
            cernctl=args.cernctl,
            cern_remote=args.cern_remote,
            eos_base=additional_eos_base,
            local_eos=additional_local_eos,
            remote_stage=f"{args.remote_stage}.additional{index}",
            include_outputs=include_outputs,
            renew_kerberos=False,
        )


def read_completed(state_path: Path) -> set[int]:
    if not state_path.exists():
        return set()
    done: set[int] = set()
    with state_path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            done.add(int(row["milestone_pct"]))
    return done


def highest_pending_milestone(
    milestones: list[int], completion: int, completed: set[int]
) -> int | None:
    completed_through = max(completed, default=-1)
    eligible = [
        milestone
        for milestone in milestones
        if completed_through < milestone <= completion
    ]
    return max(eligible, default=None)


def append_state(state_path: Path, row: dict[str, str | int]) -> None:
    state_path.parent.mkdir(parents=True, exist_ok=True)
    exists = state_path.exists()
    with state_path.open("a", newline="") as handle:
        fieldnames = [
            "date", "milestone_pct", "aa_success", "pp_success",
            "aa_events", "pp_events", "analysis_dir",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        if not exists:
            writer.writeheader()
        writer.writerow(row)


def latex_escape(value: str) -> str:
    return (
        value.replace("\\", r"\textbackslash{}")
        .replace("_", r"\_")
        .replace("%", r"\%")
        .replace("&", r"\&")
        .replace("#", r"\#")
    )


def status_block(
    *,
    campaign: str,
    milestone: int,
    aa_success: int,
    pp_success: int,
    aa_chunks: int,
    pp_chunks: int,
    aa_events_done: int,
    pp_events_done: int,
    aa_total_events: int,
    pp_total_events: int,
    updated: str,
) -> str:
    return "\n".join([
        BEGIN_MARKER,
        rf"\newcommand{{\ooAutoCampaign}}{{{latex_escape(campaign)}}}",
        rf"\newcommand{{\ooAutoMilestone}}{{{milestone}\%}}",
        rf"\newcommand{{\ooAutoAAJobs}}{{{aa_success}/{aa_chunks}}}",
        rf"\newcommand{{\ooAutoPPJobs}}{{{pp_success}/{pp_chunks}}}",
        rf"\newcommand{{\ooAutoAAEvents}}{{{aa_events_done:,}/{aa_total_events:,}}}",
        rf"\newcommand{{\ooAutoPPEvents}}{{{pp_events_done:,}/{pp_total_events:,}}}",
        rf"\newcommand{{\ooAutoUpdated}}{{{latex_escape(updated)}}}",
        rf"\newcommand{{\ooAutoFigure}}{{{DEFAULT_FIGURE}}}",
        END_MARKER,
    ])


def replace_status_block(tex_path: Path, block: str) -> None:
    text = tex_path.read_text()
    pattern = re.compile(rf"{re.escape(BEGIN_MARKER)}.*?{re.escape(END_MARKER)}", re.S)
    if not pattern.search(text):
        raise RuntimeError(f"{tex_path} does not contain OO auto-status markers")
    tex_path.write_text(pattern.sub(lambda _: block, text))


def run_analyzer(
    work: Path,
    local_eos: Path,
    additional_aa_local_eos: list[Path],
    pp_local_eos: Path,
    out_dir: Path,
    aa_events: int,
) -> None:
    analyzer = work / "analysis/analyze_oo_prehydro_pair.py"
    command = [
        "python3", str(analyzer),
        "--local-eos", str(local_eos),
        "--pp-local-eos", str(pp_local_eos),
        "--out-dir", str(out_dir),
        "--require-paired-aa",
        "--aa-events-per-chunk", str(aa_events),
    ]
    for source in additional_aa_local_eos:
        command.extend(["--additional-aa-local-eos", str(source)])
    subprocess.run(
        command,
        check=True,
    )


def read_analyzed_aa_chunks(out_dir: Path) -> int:
    source_table = out_dir / "aa_sources.tsv"
    with source_table.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"analyzer wrote no AA source rows to {source_table}")
    try:
        counts = [int(row["chunks_used"]) for row in rows]
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"invalid analyzer AA source table: {source_table}") from error
    if any(count < 0 for count in counts):
        raise ValueError(f"negative AA chunk count in {source_table}")
    return sum(counts)


def copy_figures(out_dir: Path, overleaf: Path, milestone: int) -> list[Path]:
    figure_dir = overleaf / "report/figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    copied: list[Path] = []
    prefix = "20260709-oo-c0-5-paired-prehydro-raa"
    for suffix in ("pdf", "png", "tsv"):
        source = out_dir / f"oo5360_c0_5_prehydro_overlay_raa.{suffix}"
        if not source.exists():
            continue
        current = figure_dir / f"{prefix}-current.{suffix}"
        milestone_path = figure_dir / f"{prefix}-{milestone:03d}pct.{suffix}"
        shutil.copy2(source, current)
        shutil.copy2(source, milestone_path)
        copied.extend([current, milestone_path])
    for source_name, label in (
        ("aa_sources.tsv", "sources"),
        ("aa_rejections.tsv", "rejections"),
    ):
        source = out_dir / source_name
        if not source.exists():
            continue
        current = figure_dir / f"{prefix}-current-{label}.tsv"
        milestone_path = figure_dir / f"{prefix}-{milestone:03d}pct-{label}.tsv"
        shutil.copy2(source, current)
        shutil.copy2(source, milestone_path)
        copied.extend([current, milestone_path])
    return copied


def commit_overleaf(overleaf: Path, paths: list[Path], message: str, push: bool) -> None:
    rel_paths = [str(path.relative_to(overleaf)) for path in paths if path.exists()]
    if not rel_paths:
        return
    subprocess.run(["git", "-C", str(overleaf), "add", *rel_paths], check=True)
    diff = subprocess.run(["git", "-C", str(overleaf), "diff", "--cached", "--quiet"])
    if diff.returncode == 0:
        return
    subprocess.run(["git", "-C", str(overleaf), "commit", "-m", message], check=True)
    if push:
        subprocess.run(["git", "-C", str(overleaf), "push"], check=True)


def monitor_once(args: argparse.Namespace) -> bool:
    eos_base = args.eos_base or f"/eos/user/y/yjlee/{args.campaign}"
    local_eos = args.local_eos or (args.work / "eos_snapshots" / args.campaign)
    pp_local_eos = args.pp_local_eos or local_eos
    analysis_root = args.analysis_out_root or (args.work / "analysis_1m_milestones")
    state_path = args.work / f"monitor_{args.campaign}_state.tsv"
    tex_path = args.overleaf / args.tex
    milestones = [int(item) for item in args.milestones.split(",") if item.strip()]
    aa_chunks_total = args.aa_chunks + sum(args.additional_aa_chunks)

    if args.copy:
        copy_from_eos(eos_base, local_eos)
    if args.sync_via_cernctl:
        sync_configured_aa_sources(
            args,
            eos_base=eos_base,
            local_eos=local_eos,
            include_outputs=False,
        )

    aa_status_success = len(successful_chunks(local_eos, "aa")) + sum(
        len(successful_chunks(source, "aa"))
        for source in args.additional_aa_local_eos
    )
    pp_status_success = len(successful_chunks(pp_local_eos, "pp"))
    status_completion = int(
        100.0 * min(
            aa_status_success / aa_chunks_total,
            pp_status_success / args.pp_chunks,
        )
    )
    done = read_completed(state_path)
    milestone = highest_pending_milestone(milestones, status_completion, done)
    if milestone is None:
        print(
            f"status completion={status_completion}% "
            f"aa={aa_status_success}/{aa_chunks_total} "
            f"pp={pp_status_success}/{args.pp_chunks}; no new milestone"
        )
        return False

    if args.sync_via_cernctl:
        sync_configured_aa_sources(
            args,
            eos_base=eos_base,
            local_eos=local_eos,
            include_outputs=True,
        )

    aa_success = len(available_chunks(local_eos, "aa")) + sum(
        len(available_chunks(source, "aa"))
        for source in args.additional_aa_local_eos
    )
    pp_success = len(available_chunks(pp_local_eos, "pp"))
    completion = int(
        100.0 * min(aa_success / aa_chunks_total, pp_success / args.pp_chunks)
    )
    if milestone > completion:
        print(
            f"status reached {status_completion}%, but complete local status/output "
            f"pairs are at {completion}% (aa={aa_success}/{aa_chunks_total}, "
            f"pp={pp_success}/{args.pp_chunks})"
        )
        return False

    pp_done = min(pp_success * args.pp_events, args.pp_chunks * args.pp_events)
    stamp = datetime.now().astimezone().isoformat(timespec="seconds")
    out_dir = analysis_root / f"{milestone:03d}pct"
    run_analyzer(
        args.work,
        local_eos,
        args.additional_aa_local_eos,
        pp_local_eos,
        out_dir,
        args.aa_events,
    )
    analyzed_aa_success = read_analyzed_aa_chunks(out_dir)
    analyzed_completion = int(
        100.0
        * min(
            analyzed_aa_success / aa_chunks_total,
            pp_success / args.pp_chunks,
        )
    )
    if milestone > analyzed_completion:
        print(
            f"analyzer accepted only {analyzed_aa_success}/{aa_chunks_total} "
            f"AA pairs ({analyzed_completion}%); milestone {milestone}% "
            "was not published"
        )
        return False
    aa_done = min(
        analyzed_aa_success * args.aa_events,
        aa_chunks_total * args.aa_events,
    )
    copied = copy_figures(out_dir, args.overleaf, milestone)
    block = status_block(
        campaign=args.campaign,
        milestone=milestone,
        aa_success=analyzed_aa_success,
        pp_success=pp_success,
        aa_chunks=aa_chunks_total,
        pp_chunks=args.pp_chunks,
        aa_events_done=aa_done,
        pp_events_done=pp_done,
        aa_total_events=aa_chunks_total * args.aa_events,
        pp_total_events=args.pp_chunks * args.pp_events,
        updated=stamp,
    )
    replace_status_block(tex_path, block)
    append_state(state_path, {
        "date": stamp,
        "milestone_pct": milestone,
        "aa_success": analyzed_aa_success,
        "pp_success": pp_success,
        "aa_events": aa_done,
        "pp_events": pp_done,
        "analysis_dir": str(out_dir),
    })
    if args.commit_overleaf:
        commit_overleaf(
            args.overleaf,
            [tex_path, *copied],
            f"Update OO 1M RAA milestone {milestone}pct",
            args.push_overleaf,
        )
    print(
        f"updated milestone {milestone}%: "
        f"aa={analyzed_aa_success}/{aa_chunks_total} "
        f"pp={pp_success}/{args.pp_chunks}"
    )
    return True


def main() -> int:
    args = parse_args()
    while True:
        monitor_once(args)
        if args.once:
            break
        time.sleep(args.sleep_s)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
