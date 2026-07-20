#!/usr/bin/env python3
"""Validate and summarize an MMLHI validation output directory."""

from __future__ import annotations

import argparse
import hashlib
import math
import re
from pathlib import Path


HEAVY_RE = re.compile(
    r"Heavy-quark energy-loss diagnostics:"
    r" n_heavy_steps= (?P<heavy>\d+)"
    r" n_baseline_steps= (?P<baseline>\d+)"
    r" n_drag_steps= (?P<drag>\d+)"
    r" n_diffusion_steps= (?P<diffusion>\d+)"
    r" n_diffusion_only_steps= (?P<diffusion_only>\d+)"
    r" n_invalid_steps= (?P<invalid>\d+)"
)
HARD_HEAVY_RE = re.compile(
    r"n_hard_scattered_heavy= (?P<hard_heavy>\d+)"
    r" n_hard_heavy_mass_shell_failures= (?P<hard_mass_failures>\d+)"
    r" n_hard_heavy_momentum_closure_failures= (?P<hard_closure_failures>\d+)"
    r" max_hard_heavy_mass_shell_residual= (?P<hard_max_mass_residual>[-+0-9.eE]+)"
    r" max_hard_heavy_momentum_closure_residual= (?P<hard_max_closure_residual>[-+0-9.eE]+)"
)
DYNAMIC_RE = re.compile(
    r"Dynamic unresolved Moliere diagnostics:"
    r" n_unresolved_segments_dynamic= (?P<segments>\d+)"
    r" n_unresolved_candidate_scatters= (?P<candidates>\d+)"
    r" n_unresolved_coherent_scatters= (?P<coherent>\d+)"
    r" n_unresolved_resolving_scatters= (?P<resolving>\d+)"
)
RECURSIVE_RE = re.compile(
    r"Recursive unresolved Moliere diagnostics:.*"
    r" n_recursive_frontier_permutation_checks= (?P<permutation_checks>\d+)"
    r" n_recursive_frontier_permutation_mismatches= (?P<permutation_mismatches>\d+)"
)
FATAL_PATTERNS = (
    "Both charms!?",
    "No match for charm!?",
    "Elscat problem!",
    "HYBRID configuration/runtime error:",
    "TAU Not a number",
    "Got crazy kick",
)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def output_count(path: Path, prefix: str) -> int:
    with path.open(encoding="utf-8", errors="replace") as stream:
        return sum(1 for line in stream if line.startswith(prefix))


def parse_log(path: Path) -> dict[str, int | str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    failures = [pattern for pattern in FATAL_PATTERNS if pattern in text]
    heavy_match = HEAVY_RE.search(text)
    hard_heavy_match = HARD_HEAVY_RE.search(text)
    dynamic_match = DYNAMIC_RE.search(text)
    recursive_match = RECURSIVE_RE.search(text)
    result: dict[str, int | str] = {
        "failures": ",".join(failures),
        "hadronization_retries": text.count(
            "Pythia::forceHadronLevel: hadronLevel failed; try again"
        ),
        "hadronization_giveups": text.count(
            "Pythia::forceHadronLevel: hadronLevel failed; giving up"
        ),
    }
    if heavy_match:
        result.update({key: int(value) for key, value in heavy_match.groupdict().items()})
    if hard_heavy_match:
        for key, value in hard_heavy_match.groupdict().items():
            result[key] = float(value) if "residual" in key else int(value)
    if dynamic_match:
        result.update({key: int(value) for key, value in dynamic_match.groupdict().items()})
    if recursive_match:
        result.update({key: int(value) for key, value in recursive_match.groupdict().items()})
    return result


def summarize_run(run_dir: Path) -> dict[str, int | str]:
    log = run_dir / "run.log"
    hadrons = run_dir / "result_Hadrons.out"
    partons = run_dir / "result_Partons.out"
    for required in (log, hadrons, partons):
        if not required.is_file():
            raise RuntimeError(f"missing validation output: {required}")
    result = parse_log(log)
    result.update(
        {
            "events": output_count(hadrons, "# event "),
            "hadron_end": output_count(hadrons, "end"),
            "parton_end": output_count(partons, "end"),
            "hadron_sha256": digest(hadrons),
            "parton_sha256": digest(partons),
        }
    )
    heavy_hard_scattered = 0
    with partons.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) == 6:
                numeric = [float(value) for value in fields]
                if not all(math.isfinite(value) for value in numeric):
                    raise RuntimeError(f"{run_dir.name}: non-finite parton output")
            if len(fields) == 6 and fields[4] in {"4", "-4", "5", "-5"} and fields[5] == "1":
                heavy_hard_scattered += 1
    with hadrons.open(encoding="utf-8", errors="replace") as stream:
        for line in stream:
            fields = line.split()
            if len(fields) == 6:
                numeric = [float(value) for value in fields]
                if not all(math.isfinite(value) for value in numeric):
                    raise RuntimeError(f"{run_dir.name}: non-finite hadron output")
    result["heavy_hard_scattered"] = heavy_hard_scattered
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--expected-events", type=int)
    args = parser.parse_args()

    rows: list[tuple[str, dict[str, int | str]]] = []
    for log in sorted(args.output_dir.glob("*/run.log")):
        run_dir = log.parent
        if run_dir.name.startswith("invalid_"):
            continue
        row = summarize_run(run_dir)
        if row["failures"]:
            raise RuntimeError(f"{run_dir.name}: fatal diagnostics: {row['failures']}")
        if row.get("invalid", 0) != 0:
            raise RuntimeError(f"{run_dir.name}: nonzero invalid heavy steps")
        if args.expected_events is not None and row["events"] != args.expected_events:
            raise RuntimeError(
                f"{run_dir.name}: expected {args.expected_events} events, found {row['events']}"
            )
        if row["events"] != row["hadron_end"] or row["events"] != row["parton_end"]:
            raise RuntimeError(f"{run_dir.name}: incomplete output event blocks")
        if "candidates" in row and row["candidates"] != row["coherent"] + row["resolving"]:
            raise RuntimeError(f"{run_dir.name}: unresolved-candidate accounting does not close")
        if row.get("permutation_mismatches", 0) != 0:
            raise RuntimeError(f"{run_dir.name}: Mode E permutation mismatch")
        rows.append((run_dir.name, row))

    columns = (
        "run",
        "events",
        "heavy",
        "baseline",
        "drag",
        "diffusion",
        "diffusion_only",
        "invalid",
        "hadronization_retries",
        "hadronization_giveups",
        "hard_heavy",
        "hard_mass_failures",
        "hard_closure_failures",
        "hard_max_mass_residual",
        "hard_max_closure_residual",
        "segments",
        "candidates",
        "coherent",
        "resolving",
        "permutation_checks",
        "permutation_mismatches",
        "heavy_hard_scattered",
        "hadron_sha256",
        "parton_sha256",
    )
    print("\t".join(columns))
    for name, row in rows:
        values = {"run": name, **row}
        print("\t".join(str(values.get(column, "")) for column in columns))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
