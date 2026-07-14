#!/usr/bin/env python3
"""Freeze a strict matched OO snapshot for alpha=0.37/0.335 comparison."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
import time


SUPPORT = Path(__file__).resolve().parents[1] / "cern_support"
sys.path.insert(0, str(SUPPORT))
from supervise_oo_prehydro_only import (  # noqa: E402
    load_assignments,
    ready_ids as prehydro_only_ready_ids,
)
from supervise_oo_v2 import ready_ids as paired_ready_ids  # noqa: E402


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def write_metadata_tsv(path: Path, metadata: dict[str, object]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["field", "value"])
        writer.writerows((key, value) for key, value in metadata.items())


def freeze(args: argparse.Namespace) -> dict[str, object]:
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"{output} is not empty")
    reference_sources = [path.resolve() for path in args.reference_source]
    if len(set(reference_sources)) != len(reference_sources):
        raise ValueError("reference sources must be distinct")
    alpha_source = args.alpha_source.resolve()
    manifest = args.manifest.resolve()
    assignments = load_assignments(manifest, args.target)

    reference_owner: dict[int, Path] = {}
    reference_counts: dict[str, int] = {}
    for source in reference_sources:
        ready = paired_ready_ids(source)
        overlap = set(reference_owner) & ready
        if overlap:
            preview = ",".join(str(task_id) for task_id in sorted(overlap)[:20])
            raise ValueError(f"reference sources overlap on task IDs: {preview}")
        reference_owner.update((task_id, source) for task_id in ready)
        reference_counts[str(source)] = len(ready)

    alpha_ready = prehydro_only_ready_ids(
        alpha_source,
        assignments,
        expected_alpha=args.expected_alpha,
        expected_broadening_k=args.expected_broadening_k,
    )
    selected = sorted(set(reference_owner) & alpha_ready)
    if not selected:
        raise RuntimeError("strict three-way intersection is empty")

    output.mkdir(parents=True, exist_ok=True)
    reference_out = output / "reference_eos"
    alpha_out = output / "alpha0335_eos"
    for task_id in selected:
        reference = reference_owner[task_id]
        for member, suffix in (("status", ".txt"), ("outputs", ".tar.gz")):
            link_file(
                reference / member / "aa" / f"chunk_{task_id}{suffix}",
                reference_out / member / "aa" / f"chunk_{task_id}{suffix}",
            )
            link_file(
                alpha_source / member / "aa" / f"chunk_{task_id}{suffix}",
                alpha_out / member / "aa" / f"chunk_{task_id}{suffix}",
            )

    accepted = output / "accepted_task_ids.txt"
    accepted.write_text("".join(f"{task_id}\n" for task_id in selected))
    frozen_manifest = output / "aa_task_manifest.tsv"
    shutil.copy2(manifest, frozen_manifest)
    hydro_slots = {int(assignments[task_id]["hydro_slot"]) for task_id in selected}
    metadata: dict[str, object] = {
        "snapshot_boundary": datetime.now().astimezone().isoformat(timespec="seconds"),
        "selection": "strict reference-pair AND strict alpha=0.335 completion intersection",
        "publication_status": "PROVISIONAL_COMPLETION_ORDER_MATCHED_DIAGNOSTIC",
        "target_tasks": args.target,
        "strict_reference_tasks": len(reference_owner),
        "strict_reference_source_counts": reference_counts,
        "strict_alpha0335_tasks": len(alpha_ready),
        "strict_matched_triplets": len(selected),
        "min_task_id": min(selected),
        "max_task_id": max(selected),
        "unique_hydro_slots": len(hydro_slots),
        "expected_alpha": args.expected_alpha,
        "expected_broadening_k": args.expected_broadening_k,
        "task_id_list_sha256": sha256(accepted),
        "task_manifest_sha256": sha256(frozen_manifest),
        "reference_sources": [str(path) for path in reference_sources],
        "alpha_source": str(alpha_source),
    }
    write_metadata_tsv(output / "snapshot_metadata.tsv", metadata)
    (output / "snapshot_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    )
    return metadata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-source", action="append", type=Path, required=True)
    parser.add_argument("--alpha-source", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--target", type=int, default=100_000)
    parser.add_argument("--expected-alpha", type=float, default=0.335)
    parser.add_argument("--expected-broadening-k", type=float, default=15.0)
    return parser.parse_args()


def main() -> int:
    metadata = freeze(parse_args())
    print(json.dumps(metadata, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
