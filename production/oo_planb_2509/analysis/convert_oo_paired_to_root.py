#!/usr/bin/env python3
"""Convert strict paired OO HYBRID archives to one wake-aware ROOT file."""

from __future__ import annotations

import argparse
import csv
import dataclasses
import datetime as dt
import hashlib
import io
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shlex
import struct
import subprocess
import sys
import tarfile
from typing import BinaryIO, Iterable


PAIR_HEADER = struct.Struct("<8sQiqiqiddddII")
PARTICLE_RECORD = struct.Struct("<ddddii")
PAIR_MAGIC = b"OOPAIR1\0"
CHUNK_PATTERN = re.compile(r"chunk_(\d+)\.(?:tar\.gz|txt)$")
ALLOWED_LABELS = {-2, 0, 1, 2, 3}
SCHEMA_VERSION = "oo-paired-root-v6"
JET_RADIUS_DIGITS = (1, 2, 4, 8)


class ArchiveValidationError(ValueError):
    pass


@dataclasses.dataclass(frozen=True)
class Particle:
    px: float
    py: float
    pz: float
    mass: float
    pdg_id: int
    raw_label: int


@dataclasses.dataclass(frozen=True)
class HybridEvent:
    event_number: int
    event_weight: float
    sigma_gen: float
    hard_x: float
    hard_y: float
    particles: tuple[Particle, ...]


@dataclasses.dataclass(frozen=True)
class PairedArchive:
    seed: int
    hydro_index: int
    hydro_event_id: int | None
    hydro_ncoll: int | None
    hydro_payload_sha256: str | None
    no_prehydro: HybridEvent
    with_prehydro: HybridEvent


@dataclasses.dataclass(frozen=True)
class PrehydroOnlyArchive:
    task_id: int
    seed: int
    hydro_index: int
    hydro_event_id: int
    hydro_ncoll: int
    hydro_payload_sha256: str
    energy_loss_alpha: float
    broadening_k: float
    event: HybridEvent


@dataclasses.dataclass(frozen=True)
class Source:
    name: str
    path: Path


def parse_key_value_text(text: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in text.splitlines():
        if "=" not in raw_line:
            continue
        key, value = raw_line.split("=", 1)
        values[key.strip()] = value.strip()
    return values


def parse_tsv(text: str, description: str) -> list[dict[str, str]]:
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    if reader.fieldnames is None:
        raise ArchiveValidationError(f"{description} has no TSV header")
    rows = [dict(row) for row in reader]
    if not rows:
        raise ArchiveValidationError(f"{description} has no data rows")
    return rows


def parse_hybrid_event(payload: bytes, member_name: str) -> HybridEvent:
    event_number: int | None = None
    event_weight: float | None = None
    sigma_gen: float | None = None
    hard_x: float | None = None
    hard_y: float | None = None
    particles: list[Particle] = []
    saw_end = False

    for line_number, raw_line in enumerate(payload.decode("utf-8", errors="strict").splitlines(), 1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith("# event"):
            if event_number is not None:
                raise ArchiveValidationError(f"{member_name} contains more than one event")
            fields = line.split()
            if len(fields) != 3:
                raise ArchiveValidationError(f"{member_name}:{line_number}: malformed event header")
            event_number = int(fields[2])
            continue
        if line.startswith("weight"):
            fields = line.split()
            if fields[0::2] != ["weight", "cross", "X", "Y"] or len(fields) != 8:
                raise ArchiveValidationError(f"{member_name}:{line_number}: malformed weight header")
            if event_weight is not None:
                raise ArchiveValidationError(f"{member_name} contains duplicate weight headers")
            event_weight = float(fields[1])
            sigma_gen = float(fields[3])
            hard_x = float(fields[5])
            hard_y = float(fields[7])
            continue
        if line == "end":
            saw_end = True
            continue
        if saw_end:
            raise ArchiveValidationError(f"{member_name}:{line_number}: content follows end marker")
        fields = line.split()
        if len(fields) != 6:
            raise ArchiveValidationError(f"{member_name}:{line_number}: expected six particle columns")
        px, py, pz, mass = (float(value) for value in fields[:4])
        pdg_id = int(fields[4])
        raw_label = int(fields[5])
        values = (px, py, pz, mass, event_weight, sigma_gen, hard_x, hard_y)
        if not all(value is None or math.isfinite(value) for value in values):
            raise ArchiveValidationError(f"{member_name}:{line_number}: non-finite value")
        if raw_label not in ALLOWED_LABELS:
            raise ArchiveValidationError(
                f"{member_name}:{line_number}: unsupported HYBRID label {raw_label}"
            )
        particles.append(Particle(px, py, pz, mass, pdg_id, raw_label))

    if event_number is None or event_weight is None or sigma_gen is None or hard_x is None or hard_y is None:
        raise ArchiveValidationError(f"{member_name} is missing its event or weight header")
    if not saw_end:
        raise ArchiveValidationError(f"{member_name} is missing its end marker")
    if not all(math.isfinite(value) for value in (event_weight, sigma_gen, hard_x, hard_y)):
        raise ArchiveValidationError(f"{member_name} has non-finite event metadata")
    if not particles:
        raise ArchiveValidationError(f"{member_name} has no particle records")
    return HybridEvent(event_number, event_weight, sigma_gen, hard_x, hard_y, tuple(particles))


def close_enough(first: float, second: float) -> bool:
    return math.isclose(first, second, rel_tol=1e-12, abs_tol=1e-12)


def validate_pair_metadata(no_prehydro: HybridEvent, with_prehydro: HybridEvent) -> None:
    fields = (
        ("event number", no_prehydro.event_number, with_prehydro.event_number),
        ("event weight", no_prehydro.event_weight, with_prehydro.event_weight),
        ("sigmaGen", no_prehydro.sigma_gen, with_prehydro.sigma_gen),
        ("hard X", no_prehydro.hard_x, with_prehydro.hard_x),
        ("hard Y", no_prehydro.hard_y, with_prehydro.hard_y),
    )
    for name, first, second in fields:
        equal = first == second if isinstance(first, int) else close_enough(float(first), float(second))
        if not equal:
            raise ArchiveValidationError(f"paired variants have different {name}: {first} vs {second}")


def read_tar_member(tar: tarfile.TarFile, member: tarfile.TarInfo) -> bytes:
    extracted = tar.extractfile(member)
    if extracted is None:
        raise ArchiveValidationError(f"could not extract {member.name}")
    return extracted.read()


def consistent_optional_pair_field(
    rows_by_variant: dict[str, dict[str, str]],
    key: str,
    parser: type[int] | type[str],
) -> int | str | None:
    raw_values = [
        rows_by_variant[variant].get(key, "").strip()
        for variant in sorted(rows_by_variant)
    ]
    if not any(raw_values):
        return None
    if not all(raw_values):
        raise ArchiveValidationError(f"paired variants do not both define {key}")
    try:
        values = {parser(value) for value in raw_values}
    except ValueError as error:
        raise ArchiveValidationError(f"paired variants have invalid {key}") from error
    if len(values) != 1:
        raise ArchiveValidationError(f"paired variants have different {key}")
    return values.pop()


def parse_paired_archive(path: Path, expected_chunk_id: int | None = None) -> PairedArchive:
    with tarfile.open(path, "r:gz") as tar:
        members = tar.getmembers()
        pair_summaries = [member for member in members if member.name.endswith("_pair_summary.tsv")]
        if len(pair_summaries) != 1:
            raise ArchiveValidationError(
                f"expected one pair summary, found {len(pair_summaries)}"
            )
        pair_rows = parse_tsv(
            read_tar_member(tar, pair_summaries[0]).decode("utf-8", errors="strict"),
            pair_summaries[0].name,
        )
        if len(pair_rows) != 2:
            raise ArchiveValidationError(f"pair summary has {len(pair_rows)} rows, expected two")
        rows_by_variant = {row.get("variant", ""): row for row in pair_rows}
        if set(rows_by_variant) != {"no_prehydro", "with_prehydro"}:
            raise ArchiveValidationError(f"unexpected pair-summary variants {sorted(rows_by_variant)}")

        seed_values: set[int] = set()
        hydro_values: set[int] = set()
        for variant, row in rows_by_variant.items():
            try:
                seed_values.add(int(row["seed"]))
                hydro_values.add(int(row["hydro_index"]))
                events = int(row["events"])
                return_code = int(row["returncode"])
                timeout = int(row["timeout"])
                task_id = int(row["task_id"])
            except (KeyError, TypeError, ValueError) as error:
                raise ArchiveValidationError(f"invalid {variant} pair-summary row: {error}") from error
            if events != 1 or return_code != 0 or timeout != 0:
                raise ArchiveValidationError(
                    f"{variant} summary has events={events}, returncode={return_code}, timeout={timeout}"
                )
            if expected_chunk_id is not None and task_id != expected_chunk_id:
                raise ArchiveValidationError(
                    f"pair-summary task_id={task_id} does not match chunk {expected_chunk_id}"
                )
        if len(seed_values) != 1 or len(hydro_values) != 1:
            raise ArchiveValidationError("paired variants have different seed or hydro index")
        hydro_event_id = consistent_optional_pair_field(
            rows_by_variant, "hydro_event_id", int
        )
        hydro_ncoll = consistent_optional_pair_field(
            rows_by_variant, "hydro_ncoll", int
        )
        hydro_payload_sha256 = consistent_optional_pair_field(
            rows_by_variant, "hydro_payload_sha256", str
        )

        summary_variants: dict[str, str] = {}
        for member in members:
            if not member.name.endswith("/summary.tsv"):
                continue
            rows = parse_tsv(
                read_tar_member(tar, member).decode("utf-8", errors="strict"), member.name
            )
            if len(rows) != 1 or not rows[0].get("variant"):
                raise ArchiveValidationError(f"invalid variant summary {member.name}")
            summary_variants[str(PurePosixPath(member.name).parent)] = rows[0]["variant"]

        events_by_variant: dict[str, HybridEvent] = {}
        for member in members:
            if not member.name.endswith("/HYBRID_Hadrons.out"):
                continue
            parent = str(PurePosixPath(member.name).parent)
            variant = summary_variants.get(parent)
            if variant not in {"no_prehydro", "with_prehydro"}:
                raise ArchiveValidationError(f"cannot identify variant for {member.name}")
            if variant in events_by_variant:
                raise ArchiveValidationError(f"duplicate {variant} hadron output")
            events_by_variant[variant] = parse_hybrid_event(read_tar_member(tar, member), member.name)
        if set(events_by_variant) != {"no_prehydro", "with_prehydro"}:
            raise ArchiveValidationError(
                f"archive has hadron variants {sorted(events_by_variant)}, expected both"
            )

    no_prehydro = events_by_variant["no_prehydro"]
    with_prehydro = events_by_variant["with_prehydro"]
    validate_pair_metadata(no_prehydro, with_prehydro)
    return PairedArchive(
        seed_values.pop(),
        hydro_values.pop(),
        hydro_event_id if isinstance(hydro_event_id, int) else None,
        hydro_ncoll if isinstance(hydro_ncoll, int) else None,
        hydro_payload_sha256 if isinstance(hydro_payload_sha256, str) else None,
        no_prehydro,
        with_prehydro,
    )


def parse_prehydro_only_archive(
    path: Path,
    *,
    expected_chunk_id: int | None = None,
    expected_alpha: float | None = None,
    expected_broadening_k: float | None = None,
) -> PrehydroOnlyArchive:
    with tarfile.open(path, "r:gz") as tar:
        members = tar.getmembers()
        dedicated = [
            member
            for member in members
            if member.name.endswith("_prehydro_only_summary.tsv")
        ]
        if len(dedicated) != 1:
            raise ArchiveValidationError(
                f"expected one prehydro-only summary, found {len(dedicated)}"
            )
        dedicated_rows = parse_tsv(
            read_tar_member(tar, dedicated[0]).decode("utf-8", errors="strict"),
            dedicated[0].name,
        )
        if len(dedicated_rows) != 1:
            raise ArchiveValidationError("prehydro-only summary must contain one row")

        variant_summaries = [
            member for member in members if member.name.endswith("/summary.tsv")
        ]
        if len(variant_summaries) != 1:
            raise ArchiveValidationError(
                f"expected one variant summary, found {len(variant_summaries)}"
            )
        variant_rows = parse_tsv(
            read_tar_member(tar, variant_summaries[0]).decode(
                "utf-8", errors="strict"
            ),
            variant_summaries[0].name,
        )
        if len(variant_rows) != 1:
            raise ArchiveValidationError("variant summary must contain one row")

        dedicated_row = dedicated_rows[0]
        variant_row = variant_rows[0]
        required = {
            "variant": "with_prehydro",
            "use_prehydro": "1",
            "returncode": "0",
            "timeout": "0",
        }
        for description, row in (
            ("prehydro-only summary", dedicated_row),
            ("variant summary", variant_row),
        ):
            for key, expected in required.items():
                if row.get(key) != expected:
                    raise ArchiveValidationError(
                        f"{description} has {key}={row.get(key)!r}, expected {expected!r}"
                    )

        identity_fields = (
            "task_id",
            "seed",
            "hydro_event_id",
            "hydro_ncoll",
            "hydro_payload_sha256",
            "energy_loss_alpha",
            "broadening_k",
        )
        for key in identity_fields:
            if dedicated_row.get(key) != variant_row.get(key):
                raise ArchiveValidationError(
                    f"prehydro-only summaries disagree on {key}"
                )
        dedicated_hydro_index = dedicated_row.get(
            "hydro_slot", dedicated_row.get("hydro_index")
        )
        if dedicated_hydro_index != variant_row.get("hydro_slot"):
            raise ArchiveValidationError(
                "prehydro-only summaries disagree on hydro slot/index"
            )
        if dedicated_row.get("kind") not in (None, "aa"):
            raise ArchiveValidationError("prehydro-only summary kind is not aa")
        if dedicated_row.get("events") not in (None, "1"):
            raise ArchiveValidationError("prehydro-only summary event count is not one")
        try:
            task_id = int(variant_row["task_id"])
            seed = int(variant_row["seed"])
            hydro_index = int(variant_row["hydro_slot"])
            hydro_event_id = int(variant_row["hydro_event_id"])
            hydro_ncoll = int(variant_row["hydro_ncoll"])
            payload_sha256 = variant_row["hydro_payload_sha256"]
            alpha = float(variant_row["energy_loss_alpha"])
            broadening_k = float(variant_row["broadening_k"])
        except (KeyError, TypeError, ValueError) as error:
            raise ArchiveValidationError(
                f"invalid prehydro-only summary identity: {error}"
            ) from error
        if expected_chunk_id is not None and task_id != expected_chunk_id:
            raise ArchiveValidationError(
                f"prehydro-only task_id={task_id} does not match chunk {expected_chunk_id}"
            )
        if expected_alpha is not None and not close_enough(alpha, expected_alpha):
            raise ArchiveValidationError(
                f"prehydro-only alpha={alpha}, expected {expected_alpha}"
            )
        if expected_broadening_k is not None and not close_enough(
            broadening_k, expected_broadening_k
        ):
            raise ArchiveValidationError(
                f"prehydro-only broadening K={broadening_k}, expected {expected_broadening_k}"
            )
        if not re.fullmatch(r"[0-9a-f]{64}", payload_sha256):
            raise ArchiveValidationError("invalid hydro payload SHA256")

        parent = str(PurePosixPath(variant_summaries[0].name).parent)
        hadrons = [
            member
            for member in members
            if member.name == f"{parent}/HYBRID_Hadrons.out"
        ]
        if len(hadrons) != 1 or not hadrons[0].isfile() or hadrons[0].size <= 0:
            raise ArchiveValidationError("prehydro-only archive has no unique hadron output")
        task_base = f"task_{task_id:05d}"
        forbidden = (
            f"/{task_base}/HYBRID_Hadrons.out",
            f"/{task_base}_pair_summary.tsv",
        )
        if any(
            member.name.endswith(suffix)
            for member in members
            for suffix in forbidden
        ):
            raise ArchiveValidationError("prehydro-only archive contains a baseline leg")
        event = parse_hybrid_event(read_tar_member(tar, hadrons[0]), hadrons[0].name)

    return PrehydroOnlyArchive(
        task_id=task_id,
        seed=seed,
        hydro_index=hydro_index,
        hydro_event_id=hydro_event_id,
        hydro_ncoll=hydro_ncoll,
        hydro_payload_sha256=payload_sha256,
        energy_loss_alpha=alpha,
        broadening_k=broadening_k,
        event=event,
    )


def load_aa_task_manifest(path: Path) -> dict[int, dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "task_id",
        "hard_seed",
        "hydro_slot",
        "hydro_event_id",
        "hydro_ncoll",
        "hydro_payload_sha256",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"{path}: missing required v2 task-manifest columns")
    assignments: dict[int, dict[str, str]] = {}
    for row in rows:
        task_id = int(row["task_id"])
        if task_id in assignments:
            raise ValueError(f"{path}: duplicate task_id {task_id}")
        assignments[task_id] = row
    return assignments


def validate_pair_task_assignment(
    pair: PairedArchive,
    *,
    task_id: int,
    assignment: dict[str, str],
) -> None:
    expected_values: dict[str, int | str] = {
        "seed": int(assignment["hard_seed"]),
        "hydro_index": int(assignment["hydro_slot"]),
        "hydro_event_id": int(assignment["hydro_event_id"]),
        "hydro_ncoll": int(assignment["hydro_ncoll"]),
        "hydro_payload_sha256": assignment["hydro_payload_sha256"],
    }
    for field, expected in expected_values.items():
        actual = getattr(pair, field)
        if actual != expected:
            raise ArchiveValidationError(
                f"task {task_id}: {field}={actual}, expected {expected}"
            )


def validate_prehydro_only_task_assignment(
    archive: PrehydroOnlyArchive,
    *,
    task_id: int,
    assignment: dict[str, str],
) -> None:
    expected_values: dict[str, int | str] = {
        "task_id": task_id,
        "seed": int(assignment["hard_seed"]),
        "hydro_index": int(assignment["hydro_slot"]),
        "hydro_event_id": int(assignment["hydro_event_id"]),
        "hydro_ncoll": int(assignment["hydro_ncoll"]),
        "hydro_payload_sha256": assignment["hydro_payload_sha256"],
    }
    for field, expected in expected_values.items():
        actual = getattr(archive, field)
        if actual != expected:
            raise ArchiveValidationError(
                f"task {task_id}: prehydro-only {field}={actual}, expected {expected}"
            )


def validate_hard_event_identity(
    reference: HybridEvent, candidate: HybridEvent
) -> None:
    validate_pair_metadata(reference, candidate)
    reference_markers = tuple(
        particle for particle in reference.particles if particle.raw_label == -2
    )
    candidate_markers = tuple(
        particle for particle in candidate.particles if particle.raw_label == -2
    )
    if reference_markers != candidate_markers:
        raise ArchiveValidationError(
            "matched variants have different outgoing hard-parton markers"
        )


def parse_source(value: str) -> Source:
    if "=" not in value:
        raise argparse.ArgumentTypeError("source must be NAME=LOCAL_EOS_PATH")
    name, raw_path = value.split("=", 1)
    if not name or not raw_path:
        raise argparse.ArgumentTypeError("source must have a non-empty name and path")
    path = Path(raw_path).expanduser().resolve()
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"source path is not a directory: {path}")
    return Source(name, path)


def chunk_id(path: Path) -> int:
    match = CHUNK_PATTERN.fullmatch(path.name)
    if match is None:
        raise ValueError(f"not a chunk file: {path}")
    return int(match.group(1))


def source_chunks(source: Source) -> Iterable[tuple[int, Path | None, Path | None]]:
    output_dir = source.path / "outputs" / "aa"
    status_dir = source.path / "status" / "aa"
    outputs = {chunk_id(path): path for path in output_dir.glob("chunk_*.tar.gz")}
    statuses = {chunk_id(path): path for path in status_dir.glob("chunk_*.txt")}
    for identifier in sorted(outputs.keys() | statuses.keys()):
        yield identifier, statuses.get(identifier), outputs.get(identifier)


def stable_pair_id(source_index: int, chunk: int) -> int:
    if not 0 <= source_index < 2**16:
        raise ValueError("source index does not fit pairId encoding")
    if not 0 <= chunk < 2**48:
        raise ValueError("chunk ID does not fit pairId encoding")
    return (source_index << 48) | chunk


def write_particle(stream: BinaryIO, particle: Particle) -> None:
    stream.write(
        PARTICLE_RECORD.pack(
            particle.px,
            particle.py,
            particle.pz,
            particle.mass,
            particle.pdg_id,
            particle.raw_label,
        )
    )


def write_pair(
    stream: BinaryIO,
    pair_id: int,
    source_index: int,
    chunk: int,
    pair: PairedArchive,
) -> None:
    event = pair.no_prehydro
    stream.write(
        PAIR_HEADER.pack(
            PAIR_MAGIC,
            pair_id,
            source_index,
            chunk,
            pair.hydro_index,
            pair.seed,
            event.event_number,
            event.event_weight,
            event.sigma_gen,
            event.hard_x,
            event.hard_y,
            len(pair.no_prehydro.particles),
            len(pair.with_prehydro.particles),
        )
    )
    for particle in pair.no_prehydro.particles:
        write_particle(stream, particle)
    for particle in pair.with_prehydro.particles:
        write_particle(stream, particle)


def command_output(command: list[str]) -> str:
    return subprocess.run(command, check=True, text=True, capture_output=True).stdout.strip()


def build_writer(source: Path, build_dir: Path, cxx: str, force: bool) -> Path:
    build_dir.mkdir(parents=True, exist_ok=True)
    binary = build_dir / "oo_root_tree_writer"
    if not force and binary.exists() and binary.stat().st_mtime_ns >= source.stat().st_mtime_ns:
        return binary
    root_cflags = shlex.split(command_output(["root-config", "--cflags"]))
    root_libs = shlex.split(command_output(["root-config", "--libs"]))
    fastjet_cflags = shlex.split(command_output(["fastjet-config", "--cxxflags"]))
    fastjet_libs = shlex.split(command_output(["fastjet-config", "--libs"]))
    temporary = binary.with_suffix(".tmp")
    command = [
        cxx,
        "-O3",
        "-DNDEBUG",
        "-Wall",
        "-Wextra",
        "-Wpedantic",
        *root_cflags,
        *fastjet_cflags,
        str(source),
        "-o",
        str(temporary),
        *root_libs,
        *fastjet_libs,
    ]
    subprocess.run(command, check=True)
    temporary.replace(binary)
    return binary


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(8 * 1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def validate_root(path: Path, expected_pairs: int) -> dict[str, int]:
    import uproot

    hadron_branches = {
        "pairId",
        "eventWeight",
        "hadronPt",
        "hadronEta",
        "hadronPhi",
        "hadronStatus",
        "hadronRawLabel",
        "hadronID",
    }
    jet_branches = {"pairId"}
    for radius_digit in JET_RADIUS_DIGITS:
        jet_branches.add(f"nJet{radius_digit}")
        jet_branches.update(
            f"jet{radius_digit}{suffix}"
            for suffix in (
                "Eta",
                "Phi",
                "Pt",
                "Zg",
                "Rg",
                "Mult",
                "TotalMult",
                "PtD",
                "EffectiveMultiplicity",
                "LeadingFraction",
                "NormalPtD",
                "NormalEffectiveMultiplicity",
                "LeadingNormalFraction",
                "G",
                "mass",
                "MaxKt",
                "NSD",
                "HardPartonId",
                "HardPartonPt",
                "HardPartonDR",
                "PairMatchIndex",
                "PairMatchDR",
                "PairMatchOtherPt",
                "PairMatchOtherHardPartonId",
            )
        )
        if radius_digit in (4, 8):
            jet_branches.update(
                f"jet{radius_digit}{suffix}"
                for suffix in (
                    "FormationTauF",
                    "FormationTauFSmallAngle",
                    "FormationZ",
                    "FormationTheta",
                    "FormationDeltaR",
                    "FormationKt",
                    "FormationParentE",
                    "FormationOffset",
                    "FormationInvalidSplits",
                    "FormationHardestValid",
                    "FormationHardestTauF",
                    "FormationHardestTauFSmallAngle",
                    "FormationHardestZ",
                    "FormationHardestTheta",
                    "FormationHardestDeltaR",
                    "FormationHardestKt",
                    "FormationHardestParentE",
                )
            )
    required_trees = {
        "noPrehydro/Hadrons": hadron_branches,
        "withPrehydro/Hadrons": hadron_branches,
        "noPrehydro/Jets": jet_branches,
        "withPrehydro/Jets": jet_branches,
        "Pairs": {"pairId", "sourceIndex", "chunkId", "eventWeight"},
    }
    result: dict[str, int] = {}
    with uproot.open(path) as root_file:
        for tree_name, branches in required_trees.items():
            tree = root_file[tree_name]
            if tree.num_entries != expected_pairs:
                raise RuntimeError(
                    f"{tree_name} has {tree.num_entries} entries, expected {expected_pairs}"
                )
            missing = branches - set(tree.keys())
            if missing:
                raise RuntimeError(f"{tree_name} is missing branches {sorted(missing)}")
            result[f"{tree_name}.entries"] = int(tree.num_entries)
        conversion = root_file["metadata/Conversion"].arrays(library="np")
        if int(conversion["pairCount"][0]) != expected_pairs:
            raise RuntimeError("metadata/Conversion pair count disagrees with trees")
        for name in (
            "pairCount",
            "noPrehydroHadronCount",
            "withPrehydroHadronCount",
            "noPrehydroJet1Count",
            "withPrehydroJet1Count",
            "noPrehydroJet2Count",
            "withPrehydroJet2Count",
            "noPrehydroJet4Count",
            "withPrehydroJet4Count",
            "noPrehydroJet8Count",
            "withPrehydroJet8Count",
        ):
            result[name] = int(conversion[name][0])
    return result


def git_head(path: Path) -> str | None:
    result = subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"], text=True, capture_output=True
    )
    return result.stdout.strip() if result.returncode == 0 else None


def make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        action="append",
        type=parse_source,
        required=True,
        help="repeat NAME=LOCAL_EOS_PATH for each distinct local AA source",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path)
    parser.add_argument("--summary-output", type=Path)
    parser.add_argument("--writer-log", type=Path)
    parser.add_argument(
        "--writer-source",
        type=Path,
        default=Path(__file__).with_name("oo_root_tree_writer.cc"),
    )
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--cxx", default=os.environ.get("CXX", "c++"))
    parser.add_argument("--force-build", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--limit", type=int, help="accepted-pair limit for smoke tests")
    parser.add_argument("--progress-every", type=int, default=250)
    parser.add_argument("--raw-jet-pt-min", type=float, default=1.0)
    parser.add_argument("--jet-abs-eta-max", type=float, default=5.0)
    parser.add_argument("--z-cut", type=float, default=0.1)
    parser.add_argument("--beta", type=float, default=0.0)
    parser.add_argument("--match-dr-fraction", type=float, default=0.5)
    parser.add_argument("--inventory", type=Path)
    parser.add_argument(
        "--aa-task-manifest",
        type=Path,
        help=(
            "optional combined v2 task manifest used to validate seed and hydro "
            "provenance across all sources"
        ),
    )
    parser.add_argument(
        "--allow-incomplete-aa-task-manifest",
        action="store_true",
        help=(
            "allow manifest rows absent from an explicitly provisional snapshot; "
            "present rows still receive full provenance and duplicate-ID validation"
        ),
    )
    return parser


def convert(args: argparse.Namespace) -> dict[str, object]:
    sources: list[Source] = args.source
    if len({source.name for source in sources}) != len(sources):
        raise ValueError("source names must be unique")
    if args.limit is not None and args.limit <= 0:
        raise ValueError("--limit must be positive")
    if args.progress_every <= 0:
        raise ValueError("--progress-every must be positive")
    task_manifest = (
        args.aa_task_manifest.expanduser().resolve()
        if args.aa_task_manifest is not None
        else None
    )
    task_assignments = (
        load_aa_task_manifest(task_manifest) if task_manifest is not None else None
    )
    if args.allow_incomplete_aa_task_manifest and task_assignments is None:
        raise ValueError(
            "--allow-incomplete-aa-task-manifest requires --aa-task-manifest"
        )

    output = args.output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    audit_output = (args.audit_output or output.with_suffix(".audit.tsv")).resolve()
    summary_output = (args.summary_output or output.with_suffix(".summary.json")).resolve()
    writer_log = (args.writer_log or output.with_suffix(".writer.log")).resolve()
    for path in (output, audit_output, summary_output, writer_log):
        if path.exists() and not args.overwrite:
            raise FileExistsError(f"refusing to overwrite {path}; pass --overwrite")

    writer = build_writer(
        args.writer_source.resolve(), args.build_dir.resolve(), args.cxx, args.force_build
    )
    partial_output = output.with_name(f".{output.name}.partial")
    partial_output.unlink(missing_ok=True)
    writer_command = [
        str(writer),
        "--output",
        str(partial_output),
        "--raw-jet-pt-min",
        str(args.raw_jet_pt_min),
        "--jet-abs-eta-max",
        str(args.jet_abs_eta_max),
        "--z-cut",
        str(args.z_cut),
        "--beta",
        str(args.beta),
        "--match-dr-fraction",
        str(args.match_dr_fraction),
    ]
    for source in sources:
        writer_command.extend(["--source", f"{source.name}\t{source.path}"])

    audit_fields = [
        "sourceIndex",
        "sourceName",
        "chunkId",
        "status",
        "archiveResult",
        "reason",
        "pairId",
        "seed",
        "hydroIndex",
        "hydroEventId",
        "hydroNcoll",
        "hydroPayloadSha256",
        "noPrehydroParticleRecords",
        "withPrehydroParticleRecords",
        "statusPath",
        "outputPath",
    ]
    accepted = 0
    accepted_task_ids: set[int] = set()
    rejected = 0
    skipped = 0
    scanned = 0
    source_counts: dict[str, dict[str, int]] = {
        source.name: {"accepted": 0, "rejected": 0, "skipped": 0, "scanned": 0}
        for source in sources
    }
    start = dt.datetime.now(dt.timezone.utc)
    audit_output.parent.mkdir(parents=True, exist_ok=True)
    writer_log.parent.mkdir(parents=True, exist_ok=True)

    with audit_output.open("w", newline="") as audit_stream, writer_log.open("w") as log_stream:
        audit_writer = csv.DictWriter(audit_stream, fieldnames=audit_fields, delimiter="\t")
        audit_writer.writeheader()
        process = subprocess.Popen(
            writer_command,
            stdin=subprocess.PIPE,
            stdout=log_stream,
            stderr=log_stream,
        )
        assert process.stdin is not None
        try:
            reached_limit = False
            for source_index, source in enumerate(sources):
                for chunk, status_path, archive_path in source_chunks(source):
                    if args.limit is not None and accepted >= args.limit:
                        reached_limit = True
                        break
                    scanned += 1
                    source_counts[source.name]["scanned"] += 1
                    status_values = (
                        parse_key_value_text(status_path.read_text(errors="replace"))
                        if status_path is not None
                        else {}
                    )
                    status = status_values.get("status", "missing")
                    row: dict[str, object] = {
                        "sourceIndex": source_index,
                        "sourceName": source.name,
                        "chunkId": chunk,
                        "status": status,
                        "archiveResult": "",
                        "reason": "",
                        "pairId": "",
                        "seed": "",
                        "hydroIndex": "",
                        "hydroEventId": "",
                        "hydroNcoll": "",
                        "hydroPayloadSha256": "",
                        "noPrehydroParticleRecords": "",
                        "withPrehydroParticleRecords": "",
                        "statusPath": status_path or "",
                        "outputPath": archive_path or "",
                    }
                    if status != "success":
                        row["archiveResult"] = "skipped"
                        row["reason"] = f"status={status}"
                        skipped += 1
                        source_counts[source.name]["skipped"] += 1
                        audit_writer.writerow(row)
                        continue
                    if archive_path is None:
                        row["archiveResult"] = "skipped"
                        row["reason"] = "successful status has no local archive"
                        skipped += 1
                        source_counts[source.name]["skipped"] += 1
                        audit_writer.writerow(row)
                        continue
                    try:
                        pair = parse_paired_archive(archive_path, expected_chunk_id=chunk)
                        if task_assignments is not None:
                            if chunk in accepted_task_ids:
                                raise ArchiveValidationError(
                                    f"task {chunk}: duplicate across AA sources"
                                )
                            assignment = task_assignments.get(chunk)
                            if assignment is None:
                                raise ArchiveValidationError(
                                    f"task {chunk}: missing from AA task manifest"
                                )
                            validate_pair_task_assignment(
                                pair, task_id=chunk, assignment=assignment
                            )
                        pair_id = stable_pair_id(source_index, chunk)
                        write_pair(process.stdin, pair_id, source_index, chunk, pair)
                    except (ArchiveValidationError, OSError, EOFError, tarfile.TarError, UnicodeError) as error:
                        row["archiveResult"] = "rejected"
                        row["reason"] = str(error).replace("\t", " ").replace("\n", " ")
                        rejected += 1
                        source_counts[source.name]["rejected"] += 1
                        audit_writer.writerow(row)
                        continue

                    row["archiveResult"] = "accepted"
                    row["pairId"] = pair_id
                    row["seed"] = pair.seed
                    row["hydroIndex"] = pair.hydro_index
                    row["hydroEventId"] = pair.hydro_event_id
                    row["hydroNcoll"] = pair.hydro_ncoll
                    row["hydroPayloadSha256"] = pair.hydro_payload_sha256
                    row["noPrehydroParticleRecords"] = len(pair.no_prehydro.particles)
                    row["withPrehydroParticleRecords"] = len(pair.with_prehydro.particles)
                    accepted += 1
                    accepted_task_ids.add(chunk)
                    source_counts[source.name]["accepted"] += 1
                    audit_writer.writerow(row)
                    if accepted % args.progress_every == 0:
                        audit_stream.flush()
                        elapsed = (dt.datetime.now(dt.timezone.utc) - start).total_seconds()
                        print(
                            f"accepted={accepted} rejected={rejected} skipped={skipped} "
                            f"elapsed={elapsed:.1f}s",
                            file=sys.stderr,
                            flush=True,
                        )
                if reached_limit:
                    break
            process.stdin.close()
            return_code = process.wait()
            if return_code != 0:
                raise RuntimeError(
                    f"ROOT writer exited with code {return_code}; inspect {writer_log}"
                )
            if (
                task_assignments is not None
                and args.limit is None
                and not args.allow_incomplete_aa_task_manifest
            ):
                missing_task_ids = set(task_assignments) - accepted_task_ids
                if missing_task_ids:
                    preview = ",".join(
                        str(value) for value in sorted(missing_task_ids)[:20]
                    )
                    raise RuntimeError(
                        "AA task-manifest closure failed: "
                        f"{len(missing_task_ids)} task(s) were not accepted; "
                        f"first IDs: {preview}"
                    )
        except BaseException:
            try:
                process.stdin.close()
            except (BrokenPipeError, OSError):
                pass
            if process.poll() is None:
                process.terminate()
            process.wait()
            partial_output.unlink(missing_ok=True)
            raise

    if accepted == 0:
        partial_output.unlink(missing_ok=True)
        raise RuntimeError("no strict paired archives were accepted")
    root_totals = validate_root(partial_output, accepted)
    partial_output.replace(output)
    finish = dt.datetime.now(dt.timezone.utc)
    source_repository = args.writer_source.resolve().parents[3]
    summary: dict[str, object] = {
        "schemaVersion": SCHEMA_VERSION,
        "startedUtc": start.isoformat(),
        "finishedUtc": finish.isoformat(),
        "elapsedSeconds": (finish - start).total_seconds(),
        "sourceRepositoryHead": git_head(source_repository),
        "sources": [
            {"sourceIndex": index, "name": source.name, "path": str(source.path)}
            for index, source in enumerate(sources)
        ],
        "sourceCounts": source_counts,
        "scannedChunkRecords": scanned,
        "acceptedPairs": accepted,
        "rejectedArchives": rejected,
        "skippedChunkRecords": skipped,
        "acceptedPairLimit": args.limit,
        "rootTotals": root_totals,
        "jetConfiguration": {
            "algorithm": "anti-kt E-scheme",
            "radii": [0.1, 0.2, 0.4, 0.8],
            "rawJetPtMinGeV": args.raw_jet_pt_min,
            "jetAbsEtaMax": args.jet_abs_eta_max,
            "constituentPtMinGeV": 0.0,
            "constituentAbsEtaMax": None,
            "softDropZCut": args.z_cut,
            "softDropBeta": args.beta,
            "pairMatchDRFraction": args.match_dr_fraction,
            "hardPartonMatchDRFraction": 1.0,
            "totalMultiplicity": (
                "NNormal + NPositiveWake - NNegativeWake, where negative wake "
                "includes raw-label-2 negative thermal hadrons and raw-label-3 "
                "hadronized holes"
            ),
            "formationTime": {
                "storedRadii": [0.4, 0.8],
                "correctedJetPtMinExclusiveGeV": 30.0,
                "reclustering": "Cambridge-Aachen E-scheme",
                "constituents": "normal plus positive wake",
                "negativeCorrection": "ghost association and jet-level 4MomSub only",
                "hbarCGeVFm": 0.19732698,
                "singleSplitting": "global maximum kT over the full C/A tree",
            },
        },
        "output": str(output),
        "outputBytes": output.stat().st_size,
        "outputSha256": sha256(output),
        "audit": str(audit_output),
        "writerLog": str(writer_log),
        "writerBinary": str(writer),
        "writerSource": str(args.writer_source.resolve()),
        "writerSourceSha256": sha256(args.writer_source.resolve()),
        "rootVersion": command_output(["root-config", "--version"]),
        "fastjetVersion": command_output(["fastjet-config", "--version"]),
        "inventory": str(args.inventory.resolve()) if args.inventory else None,
        "aaTaskManifest": (
            {
                "path": str(task_manifest),
                "sha256": sha256(task_manifest),
                "rows": len(task_assignments),
                "acceptedRows": len(accepted_task_ids),
                "missingRows": len(set(task_assignments) - accepted_task_ids),
                "closure": (
                    "PASS"
                    if args.limit is None and set(task_assignments) == accepted_task_ids
                    else (
                        "INCOMPLETE_ALLOWED"
                        if args.allow_incomplete_aa_task_manifest
                        else "NOT_REQUIRED_LIMITED"
                    )
                ),
            }
            if task_manifest is not None and task_assignments is not None
            else None
        ),
        "command": shlex.join(sys.argv),
    }
    summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> int:
    args = make_parser().parse_args()
    summary = convert(args)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
