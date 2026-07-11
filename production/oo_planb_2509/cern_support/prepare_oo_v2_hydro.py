#!/usr/bin/env python3
"""Prepare checksum-pinned OO hydro payloads and a balanced v2 task manifest."""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import random
import re
import struct
import tarfile


ZENODO_C0_5_MD5 = "32f45a6a7e7fd513b80b97857a7a136b"
EVENT_PATTERN = re.compile(r"^C0-5/hydro_results_(\d+)/([^/]+)$")


def digest_bytes(payload: bytes, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    digest.update(payload)
    return digest.hexdigest()


def digest_file(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as stream:
        while block := stream.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def parse_hydro_header(payload: bytes) -> dict[str, float | int]:
    if len(payload) < 64:
        raise ValueError("hydro evolution is too small to contain its header")
    header = struct.unpack("<16f", payload[:64])
    return {
        "tau0_fm": float(header[0]),
        "dtau_fm": float(header[1]),
        "nx": int(header[2]),
        "dx_fm": float(header[3]),
        "ny": int(header[5]),
        "dy_fm": float(header[6]),
        "n_fields": int(header[15]),
    }


def validate_hydro_header(event_id: int, header: dict[str, float | int]) -> None:
    expected_floats = {
        "tau0_fm": 0.4,
        "dtau_fm": 0.1,
        "dx_fm": 0.2,
        "dy_fm": 0.2,
    }
    expected_ints = {"nx": 150, "ny": 150, "n_fields": 11}
    for key, expected in expected_floats.items():
        if not math.isclose(float(header[key]), expected, rel_tol=0, abs_tol=1e-6):
            raise ValueError(
                f"event {event_id}: unexpected hydro {key}={header[key]}"
            )
    for key, expected in expected_ints.items():
        if int(header[key]) != expected:
            raise ValueError(
                f"event {event_id}: unexpected hydro {key}={header[key]}"
            )


def ncoll_count(payload: bytes) -> int:
    lines = payload.decode("utf-8", errors="strict").splitlines()
    count = sum(1 for line in lines if line.strip() and not line.lstrip().startswith("#"))
    if count <= 0:
        raise ValueError("Ncoll list contains no binary-collision positions")
    return count


def add_tar_bytes(archive: tarfile.TarFile, name: str, payload: bytes) -> None:
    member = tarfile.TarInfo(name)
    member.size = len(payload)
    member.mode = 0o644
    member.mtime = 0
    member.uid = 0
    member.gid = 0
    member.uname = ""
    member.gname = ""
    archive.addfile(member, io.BytesIO(payload))


def write_event_payload(path: Path, files: dict[str, bytes]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=1
        ) as compressed:
            with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as archive:
                for name in sorted(files):
                    add_tar_bytes(archive, name, files[name])


def source_event_files(archive: tarfile.TarFile) -> dict[int, set[str]]:
    events: dict[int, set[str]] = {}
    for member in archive.getmembers():
        if not member.isfile():
            continue
        match = EVENT_PATTERN.fullmatch(member.name)
        if not match:
            continue
        event_id = int(match.group(1))
        basename = match.group(2)
        event = events.setdefault(event_id, set())
        if basename in event:
            raise ValueError(f"event {event_id}: duplicate archive member {basename}")
        event.add(basename)
    return events


def read_member(archive: tarfile.TarFile, member: tarfile.TarInfo) -> bytes:
    stream = archive.extractfile(member)
    if stream is None:
        raise ValueError(f"could not extract {member.name}")
    return stream.read()


def largest_remainder_quotas(
    ncoll_values: list[int], total: int
) -> list[int]:
    if total <= 0:
        raise ValueError("allocation total must be positive")
    denominator = sum(ncoll_values)
    if denominator <= 0:
        raise ValueError("sum of Ncoll values must be positive")
    raw = [Fraction(total * value, denominator) for value in ncoll_values]
    quotas = [value.numerator // value.denominator for value in raw]
    remaining = total - sum(quotas)
    order = sorted(
        range(len(raw)),
        key=lambda index: (raw[index] - quotas[index], -index),
        reverse=True,
    )
    for index in order[:remaining]:
        quotas[index] += 1
    if sum(quotas) != total:
        raise AssertionError("largest-remainder allocation does not close")
    return quotas


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def make_task_rows(
    hydro_rows: list[dict[str, object]],
    *,
    hard_events: int,
    block_size: int,
    seed_offset: int,
    shuffle_seed: int,
) -> tuple[list[dict[str, object]], list[int]]:
    if hard_events % block_size != 0:
        raise ValueError("hard-event total must be divisible by milestone block size")
    ncoll_values = [int(row["ncoll"]) for row in hydro_rows]
    block_quotas = largest_remainder_quotas(ncoll_values, block_size)
    if any(quota <= 0 for quota in block_quotas):
        raise ValueError("each hydro event must receive at least one task per block")
    task_rows: list[dict[str, object]] = []
    blocks = hard_events // block_size
    for block in range(blocks):
        hydro_slots = [
            slot
            for slot, quota in enumerate(block_quotas)
            for _ in range(quota)
        ]
        random.Random(shuffle_seed + block).shuffle(hydro_slots)
        for offset, slot in enumerate(hydro_slots):
            task_id = block * block_size + offset
            hydro = hydro_rows[slot]
            task_rows.append(
                {
                    "task_id": task_id,
                    "hard_seed": seed_offset + task_id,
                    "milestone_block": block + 1,
                    "hydro_slot": slot,
                    "hydro_event_id": hydro["event_id"],
                    "hydro_ncoll": hydro["ncoll"],
                    "hydro_dir": hydro["hydro_dir"],
                    "hydro_payload_key": hydro["payload_key"],
                    "hydro_payload_sha256": hydro["payload_sha256"],
                }
            )
    if len(task_rows) != hard_events:
        raise AssertionError("task manifest has the wrong number of rows")
    return task_rows, block_quotas


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--expected-events", type=int, default=500)
    parser.add_argument("--expected-md5", default=ZENODO_C0_5_MD5)
    parser.add_argument("--hard-events", type=int, default=50_000)
    parser.add_argument("--milestone-block-size", type=int, default=5_000)
    parser.add_argument("--seed-offset", type=int, default=900_000)
    parser.add_argument("--shuffle-seed", type=int, default=20260711)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    archive_path = args.archive.resolve()
    out_dir = args.out_dir.resolve()
    if not archive_path.is_file():
        raise FileNotFoundError(archive_path)
    actual_md5 = digest_file(archive_path, "md5")
    if actual_md5 != args.expected_md5:
        raise ValueError(
            f"{archive_path}: expected MD5 {args.expected_md5}, got {actual_md5}"
        )
    if out_dir.exists() and any(out_dir.iterdir()) and not args.overwrite:
        raise FileExistsError(f"{out_dir} is not empty; pass --overwrite")
    out_dir.mkdir(parents=True, exist_ok=True)
    payload_dir = out_dir / "hydro_payloads"
    payload_dir.mkdir(parents=True, exist_ok=True)

    with tarfile.open(archive_path, "r:gz") as source:
        event_files = source_event_files(source)
    if len(event_files) != args.expected_events:
        raise ValueError(
            f"expected {args.expected_events} hydro events, found {len(event_files)}"
        )
    event_ids = sorted(event_files)
    slot_by_event = {event_id: slot for slot, event_id in enumerate(event_ids)}
    for event_id, basenames in event_files.items():
        required = {
            "evolution_all_xyeta.dat",
            f"NcollList{event_id}.dat",
            "music_input",
            "run.log",
        }
        missing = required - basenames
        if missing:
            raise ValueError(f"event {event_id}: missing {sorted(missing)}")

    hydro_rows_by_event: dict[int, dict[str, object]] = {}

    def finalize_event(event_id: int, files: dict[str, bytes]) -> None:
        if event_id in hydro_rows_by_event:
            raise ValueError(f"event {event_id}: members are not contiguous in archive")
        slot = slot_by_event[event_id]
        evolution_name = "evolution_all_xyeta.dat"
        ncoll_name = f"NcollList{event_id}.dat"
        required = (evolution_name, ncoll_name, "music_input", "run.log")
        for name in required:
            if name not in files:
                raise ValueError(f"event {event_id}: streaming pass missed {name}")
        evolution = files[evolution_name]
        ncoll_payload = files[ncoll_name]
        music_input = files["music_input"]
        run_log = files["run.log"]
        ncoll = ncoll_count(ncoll_payload)
        header = parse_hydro_header(evolution)
        validate_hydro_header(event_id, header)
        hydro_dir = f"C0-5_event_{event_id:05d}"
        payload_name = f"event_{event_id:05d}.tar.gz"
        payload_path = payload_dir / payload_name
        evolution_sha256 = digest_bytes(evolution)
        ncoll_sha256 = digest_bytes(ncoll_payload)
        readme = (
            "# OO 5.36 TeV v2 staged hydro event\n\n"
            "centrality = C0-5\n"
            f"hydro_slot = {slot}\n"
            f"event_id = {event_id}\n"
            f"ncoll_positions = {ncoll}\n"
            f"source_archive_md5 = {actual_md5}\n"
            f"source_hydro = C0-5/hydro_results_{event_id}/evolution_all_xyeta.dat\n"
            f"source_ncoll = C0-5/hydro_results_{event_id}/{ncoll_name}\n"
            f"evolution_sha256 = {evolution_sha256}\n"
            f"ncoll_sha256 = {ncoll_sha256}\n"
            f"tau0_fm = {header['tau0_fm']:.9g}\n"
            f"dtau_fm = {header['dtau_fm']:.9g}\n"
            f"nx = {header['nx']}\n"
            f"ny = {header['ny']}\n"
            f"dx_fm = {header['dx_fm']:.9g}\n"
            f"dy_fm = {header['dy_fm']:.9g}\n"
            f"n_fields = {header['n_fields']}\n"
        ).encode()
        write_event_payload(
            payload_path,
            {
                f"{hydro_dir}/NcollList.dat": ncoll_payload,
                f"{hydro_dir}/README_staged_event.txt": readme,
                f"{hydro_dir}/evolution_all_xyeta.dat": evolution,
                f"{hydro_dir}/music_input": music_input,
                f"{hydro_dir}/run.log": run_log,
            },
        )
        hydro_rows_by_event[event_id] = {
            "hydro_slot": slot,
            "centrality": "C0-5",
            "event_id": event_id,
            "ncoll": ncoll,
            "hydro_dir": hydro_dir,
            "payload_key": f"hydro/C0-5/{payload_name}",
            "payload_size": payload_path.stat().st_size,
            "payload_sha256": digest_file(payload_path),
            "evolution_size": len(evolution),
            "evolution_sha256": evolution_sha256,
            "ncoll_sha256": ncoll_sha256,
            "tau0_fm": f"{float(header['tau0_fm']):.9g}",
            "dtau_fm": f"{float(header['dtau_fm']):.9g}",
            "nx": header["nx"],
            "ny": header["ny"],
            "dx_fm": f"{float(header['dx_fm']):.9g}",
            "dy_fm": f"{float(header['dy_fm']):.9g}",
            "n_fields": header["n_fields"],
        }

    current_event_id: int | None = None
    current_files: dict[str, bytes] = {}
    with tarfile.open(archive_path, "r|gz") as source:
        for member in source:
            if not member.isfile():
                continue
            match = EVENT_PATTERN.fullmatch(member.name)
            if not match:
                continue
            event_id = int(match.group(1))
            basename = match.group(2)
            if current_event_id is None:
                current_event_id = event_id
            elif event_id != current_event_id:
                finalize_event(current_event_id, current_files)
                current_event_id = event_id
                current_files = {}
            if basename in {
                "evolution_all_xyeta.dat",
                f"NcollList{event_id}.dat",
                "music_input",
                "run.log",
            }:
                current_files[basename] = read_member(source, member)
    if current_event_id is not None:
        finalize_event(current_event_id, current_files)
    if set(hydro_rows_by_event) != set(event_ids):
        missing = sorted(set(event_ids) - set(hydro_rows_by_event))
        raise ValueError(f"streaming pass did not prepare all events; missing {missing}")
    hydro_rows = [hydro_rows_by_event[event_id] for event_id in event_ids]

    hydro_fields = list(hydro_rows[0])
    hydro_manifest = out_dir / "hydro_manifest.tsv"
    write_tsv(hydro_manifest, hydro_fields, hydro_rows)
    task_rows, block_quotas = make_task_rows(
        hydro_rows,
        hard_events=args.hard_events,
        block_size=args.milestone_block_size,
        seed_offset=args.seed_offset,
        shuffle_seed=args.shuffle_seed,
    )
    task_manifest = out_dir / "aa_task_manifest.tsv"
    write_tsv(task_manifest, list(task_rows[0]), task_rows)

    totals = [0] * len(hydro_rows)
    for row in task_rows:
        totals[int(row["hydro_slot"])] += 1
    summary = {
        "archive": str(archive_path),
        "archiveBytes": archive_path.stat().st_size,
        "archiveMd5": actual_md5,
        "centrality": "C0-5",
        "hydroEvents": len(hydro_rows),
        "sumNcoll": sum(int(row["ncoll"]) for row in hydro_rows),
        "minNcoll": min(int(row["ncoll"]) for row in hydro_rows),
        "maxNcoll": max(int(row["ncoll"]) for row in hydro_rows),
        "hardEvents": len(task_rows),
        "seedOffset": args.seed_offset,
        "milestoneBlockSize": args.milestone_block_size,
        "milestoneBlocks": args.hard_events // args.milestone_block_size,
        "shuffleSeed": args.shuffle_seed,
        "minTasksPerHydro": min(totals),
        "maxTasksPerHydro": max(totals),
        "payloadBytes": sum(int(row["payload_size"]) for row in hydro_rows),
        "hydroManifest": hydro_manifest.name,
        "hydroManifestSha256": digest_file(hydro_manifest),
        "taskManifest": task_manifest.name,
        "taskManifestSha256": digest_file(task_manifest),
        "blockQuotaMin": min(block_quotas),
        "blockQuotaMax": max(block_quotas),
        "allocation": "largest remainder proportional to Ncoll in every milestone block",
    }
    summary_path = out_dir / "hydro_preparation_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
