#!/usr/bin/env python3
"""Freeze and select a strict completion-order OO V2 analysis snapshot."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


CERN_SUPPORT = Path(__file__).resolve().parents[1] / "cern_support"
sys.path.insert(0, str(CERN_SUPPORT))
from supervise_oo_v2 import ready_ids  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_head(path: Path) -> str | None:
    result = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        text=True,
        capture_output=True,
        check=False,
    )
    return result.stdout.strip() if result.returncode == 0 else None


def link_file(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    for attempt in range(3):
        try:
            os.link(source, destination)
            return
        except FileNotFoundError:
            if attempt == 2:
                raise
            time.sleep(0.05)


def link_glob(source_dir: Path, destination_dir: Path, pattern: str) -> int:
    count = 0
    for source in sorted(source_dir.glob(pattern)):
        link_file(source, destination_dir / source.name)
        count += 1
    return count


def load_manifest(path: Path) -> tuple[list[str], dict[int, dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    required = {"task_id", "hydro_slot"}
    if not rows or not required.issubset(fieldnames):
        raise ValueError(f"{path}: missing required fields {sorted(required)}")
    assignments: dict[int, dict[str, str]] = {}
    for row in rows:
        task_id = int(row["task_id"])
        if task_id in assignments:
            raise ValueError(f"{path}: duplicate task_id {task_id}")
        assignments[task_id] = row
    return fieldnames, assignments


def write_metadata_tsv(path: Path, metadata: dict[str, object]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["field", "value"])
        writer.writerows((key, value) for key, value in metadata.items())


def freeze_snapshot(
    *,
    source: Path,
    manifest: Path,
    output: Path,
    source_repository: Path,
) -> dict[str, object]:
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"{output} is not empty")
    source = source.resolve()
    manifest = manifest.resolve()
    output.mkdir(parents=True, exist_ok=True)
    boundary = datetime.now().astimezone().isoformat(timespec="seconds")

    raw = output / "eos_snapshot_raw"
    raw_status = raw / "status/aa"
    raw_outputs = raw / "outputs/aa"
    status_count = link_glob(source / "status/aa", raw_status, "chunk_*.txt")
    output_count = link_glob(
        source / "outputs/aa", raw_outputs, "chunk_*.tar.gz"
    )

    _, assignments = load_manifest(manifest)
    selected = sorted(ready_ids(raw))
    unexpected = set(selected) - set(assignments)
    if unexpected:
        preview = ",".join(str(task_id) for task_id in sorted(unexpected)[:20])
        raise ValueError(f"strict snapshot has tasks absent from manifest: {preview}")
    if not selected:
        raise RuntimeError("frozen snapshot has no strict paired tasks")

    selected_root = output / "local_eos"
    for task_id in selected:
        link_file(
            raw_status / f"chunk_{task_id}.txt",
            selected_root / f"status/aa/chunk_{task_id}.txt",
        )
        link_file(
            raw_outputs / f"chunk_{task_id}.tar.gz",
            selected_root / f"outputs/aa/chunk_{task_id}.tar.gz",
        )

    accepted_path = output / "accepted_task_ids.txt"
    accepted_path.write_text("".join(f"{task_id}\n" for task_id in selected))
    frozen_manifest = output / "aa_task_manifest.tsv"
    shutil.copy2(manifest, frozen_manifest)
    hydro_slots = {int(assignments[task_id]["hydro_slot"]) for task_id in selected}
    metadata: dict[str, object] = {
        "snapshot_boundary": boundary,
        "strict_accepted_pairs": len(selected),
        "min_task_id": min(selected),
        "max_task_id": max(selected),
        "unique_hydro_slots": len(hydro_slots),
        "raw_status_files": status_count,
        "raw_output_files": output_count,
        "task_id_list_sha256": sha256(accepted_path),
        "task_manifest_sha256": sha256(frozen_manifest),
        "source_commit": git_head(source_repository.resolve()) or "unknown",
        "selection": "strict-complete-at-snapshot (completion-order selected)",
        "publication_status": "PROVISIONAL_DIAGNOSTIC_NOT_UNBIASED",
    }
    write_metadata_tsv(output / "snapshot_metadata.tsv", metadata)
    (output / "snapshot_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )
    return metadata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--source-repository", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metadata = freeze_snapshot(
        source=args.source,
        manifest=args.manifest,
        output=args.output,
        source_repository=args.source_repository,
    )
    print(json.dumps(metadata, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
