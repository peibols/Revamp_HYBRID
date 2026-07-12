#!/usr/bin/env python3
"""Create a nonoverlapping OO V2 continuation task manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


INTEGER_FIELDS = (
    "task_id",
    "hard_seed",
    "milestone_block",
    "hydro_slot",
    "hydro_event_id",
    "hydro_ncoll",
)


def digest_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    missing = set(INTEGER_FIELDS) - set(fieldnames)
    if not rows or missing:
        raise ValueError(f"{path}: empty manifest or missing fields {sorted(missing)}")
    for row in rows:
        for field in INTEGER_FIELDS:
            try:
                int(row[field])
            except (TypeError, ValueError) as error:
                raise ValueError(f"{path}: invalid integer field {field}") from error
    task_ids = [int(row["task_id"]) for row in rows]
    if task_ids != list(range(task_ids[0], task_ids[0] + len(task_ids))):
        raise ValueError(f"{path}: task IDs must be ordered and contiguous")
    hard_seeds = [int(row["hard_seed"]) for row in rows]
    if hard_seeds != list(range(hard_seeds[0], hard_seeds[0] + len(hard_seeds))):
        raise ValueError(f"{path}: hard seeds must be ordered and contiguous")
    return fieldnames, rows


def shifted_rows(
    rows: list[dict[str, str]],
    *,
    task_id_offset: int,
    milestone_block_offset: int,
) -> list[dict[str, str]]:
    if task_id_offset <= 0:
        raise ValueError("task ID offset must be positive")
    if milestone_block_offset <= 0:
        raise ValueError("milestone block offset must be positive")
    shifted: list[dict[str, str]] = []
    for source in rows:
        row = dict(source)
        row["task_id"] = str(int(source["task_id"]) + task_id_offset)
        row["hard_seed"] = str(int(source["hard_seed"]) + task_id_offset)
        row["milestone_block"] = str(
            int(source["milestone_block"]) + milestone_block_offset
        )
        shifted.append(row)
    return shifted


def write_manifest(
    path: Path, fieldnames: list[str], rows: list[dict[str, str]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--combined-output", type=Path)
    parser.add_argument("--summary", type=Path)
    parser.add_argument("--task-id-offset", type=int, required=True)
    parser.add_argument("--milestone-block-offset", type=int, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    fieldnames, source_rows = read_manifest(args.source)
    continuation = shifted_rows(
        source_rows,
        task_id_offset=args.task_id_offset,
        milestone_block_offset=args.milestone_block_offset,
    )
    source_ids = {int(row["task_id"]) for row in source_rows}
    continuation_ids = {int(row["task_id"]) for row in continuation}
    if source_ids & continuation_ids:
        raise ValueError("source and continuation task IDs overlap")
    write_manifest(args.output, fieldnames, continuation)
    if args.combined_output is not None:
        write_manifest(
            args.combined_output,
            fieldnames,
            [*source_rows, *continuation],
        )

    summary = {
        "source": str(args.source.resolve()),
        "sourceRows": len(source_rows),
        "sourceSha256": digest_file(args.source),
        "taskIdOffset": args.task_id_offset,
        "milestoneBlockOffset": args.milestone_block_offset,
        "continuation": str(args.output.resolve()),
        "continuationRows": len(continuation),
        "continuationTaskMin": min(continuation_ids),
        "continuationTaskMax": max(continuation_ids),
        "continuationSeedMin": int(continuation[0]["hard_seed"]),
        "continuationSeedMax": int(continuation[-1]["hard_seed"]),
        "continuationSha256": digest_file(args.output),
        "hydroAssignment": "exact source-manifest repetition with disjoint hard seeds",
    }
    if args.combined_output is not None:
        summary.update(
            {
                "combined": str(args.combined_output.resolve()),
                "combinedRows": len(source_rows) + len(continuation),
                "combinedSha256": digest_file(args.combined_output),
            }
        )
    rendered = json.dumps(summary, indent=2, sort_keys=True) + "\n"
    if args.summary is not None:
        args.summary.parent.mkdir(parents=True, exist_ok=True)
        args.summary.write_text(rendered)
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
