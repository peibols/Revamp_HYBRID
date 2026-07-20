#!/usr/bin/env python3
"""Check one recoil/hole pair for each accepted unresolved scattering."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("C", "D", "E"))
    parser.add_argument("history", type=Path)
    args = parser.parse_args()

    selected: Counter[tuple[str, str]] = Counter()
    recoilers: Counter[tuple[str, str]] = Counter()
    holes: Counter[tuple[str, str]] = Counter()
    with args.history.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            key = (row["event_id"], row["time"])
            if args.mode in {"C", "D"} and row["record_type"] == "dynamic_resolution_test":
                selected[key] += 1
            if (
                args.mode == "E"
                and row["record_type"] == "moliere_kick"
                and row["note"].startswith("modeE_")
            ):
                selected[key] += 1
            if row["record_type"] == "medium_response" and row["note"] == "recoiler":
                recoilers[key] += 1
            if row["record_type"] == "medium_response" and row["note"] == "hole":
                holes[key] += 1

    if not selected:
        raise RuntimeError(f"Mode {args.mode}: no accepted unresolved scatterings in {args.history}")
    for key, count in selected.items():
        if recoilers[key] != count or holes[key] != count:
            raise RuntimeError(
                f"Mode {args.mode} event/time {key}: {count} accepted scatters, "
                f"{recoilers[key]} recoilers, {holes[key]} holes"
            )
    print(
        f"Mode {args.mode}: {sum(selected.values())} accepted unresolved scatterings "
        "each have exactly one recoil and one hole"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
