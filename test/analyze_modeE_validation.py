#!/usr/bin/env python3
"""Summarize Mode-E recursive unresolved-Moliere validation histories.

The input is the TSV produced by dump_hybrid_evolution_history=true.  The
script is intentionally lightweight so it can be run on scratch samples without
ROOT.  It reports the concrete Dani/Krishna-style nested tests, angular changes
with respect to each shower root, and the Mode-E parent-opening closure records.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


def as_int(value: str, default: int = -1) -> int:
    try:
        return int(value)
    except Exception:
        return default


def as_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value)
    except Exception:
        return default


def dphi(a: Tuple[float, float, float, float], b: Tuple[float, float, float, float]) -> float:
    da = math.atan2(a[1], a[0]) - math.atan2(b[1], b[0])
    while da > math.pi:
        da -= 2.0 * math.pi
    while da < -math.pi:
        da += 2.0 * math.pi
    return abs(da)


def momentum(row: Dict[str, str]) -> Tuple[float, float, float, float]:
    return (as_float(row["px"]), as_float(row["py"]), as_float(row["pz"]), as_float(row["E"]))


def load(path: Path) -> List[Dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def build_parent_maps(rows: Iterable[Dict[str, str]]):
    parent: Dict[Tuple[int, int], int] = {}
    split_p: Dict[Tuple[int, int], Tuple[float, float, float, float]] = {}
    for row in rows:
        if row["record_type"] != "split":
            continue
        event = as_int(row["event_id"])
        parton = as_int(row["parton_id"])
        parent[(event, parton)] = as_int(row["parent_id"])
        split_p[(event, parton)] = momentum(row)
    return parent, split_p


def root_of(parent: Dict[Tuple[int, int], int], event: int, parton: int) -> int:
    seen = set()
    cur = parton
    while (event, cur) in parent and parent[(event, cur)] >= 0 and cur not in seen:
        seen.add(cur)
        cur = parent[(event, cur)]
    return cur


def extract_note_float(note: str, key: str) -> float | None:
    m = re.search(rf"{re.escape(key)}=([-+0-9.eE]+)", note)
    if not m:
        return None
    return as_float(m.group(1))


def summarize_case(name: str, rows: List[Dict[str, str]], max_examples: int) -> str:
    parent, split_p = build_parent_maps(rows)
    out: List[str] = []
    out.append(f"## {name}")
    out.append("")

    tests = [r for r in rows if r["record_type"] == "recursive_resolution_test"]
    groups: Dict[Tuple[int, str, int], List[Dict[str, str]]] = defaultdict(list)
    for row in tests:
        groups[(as_int(row["event_id"]), row["time"], as_int(row["parton_id"]))].append(row)
    nested = [g for g in groups.values() if len(g) >= 2]
    resolving = [g for g in groups.values() if any("resolves" in r["note"] for r in g)]
    nested_resolving = [g for g in nested if any("resolves" in r["note"] for r in g)]

    out.append("### Recursive candidate tests")
    out.append(f"- resolution-test records: {len(tests)}")
    out.append(f"- candidate groups: {len(groups)}")
    out.append(f"- nested groups with two or more bottom-up tests: {len(nested)}")
    out.append(f"- groups with a resolving test: {len(resolving)}")
    out.append(f"- nested groups with a resolving test: {len(nested_resolving)}")
    out.append("")

    examples = nested_resolving or nested
    if examples:
        out.append("### Concrete nested examples")
        for g in examples[:max_examples]:
            g = sorted(g, key=lambda r: (as_int(r["parent_id"]), as_int(r["d1"]), as_int(r["d2"])))
            first = g[0]
            event = as_int(first["event_id"])
            probe = as_int(first["parton_id"])
            time = as_float(first["time"])
            chain = []
            for r in g:
                outcome = "resolves" if "resolves" in r["note"] else "coherent"
                chain.append(
                    f"parent {r['parent_id']} tests child/sibling {r['d1']}/{r['d2']} "
                    f"qd={as_float(r['qperp']):.4g} -> {outcome}"
                )
            root = root_of(parent, event, probe)
            vac_angle = None
            if (event, probe) in split_p and (event, root) in split_p:
                vac_angle = dphi(split_p[(event, probe)], split_p[(event, root)])
            angle_text = f", vacuum DeltaPhi(probe,root)={vac_angle:.4g}" if vac_angle is not None else ""
            out.append(f"- event {event}, probe {probe}, root {root}, t={time:.6g}{angle_text}")
            for item in chain:
                out.append(f"  - {item}")
        out.append("")

    kicks = [r for r in rows if r["record_type"] == "moliere_kick" and r["note"].startswith("modeE_recursive")]
    coherent_kicks = [r for r in kicks if "coherent" in r["note"]]
    resolving_kicks = [r for r in kicks if "resolving" in r["note"]]
    angle_deltas = []
    for r in kicks:
        event = as_int(r["event_id"])
        parton = as_int(r["parton_id"])
        root = root_of(parent, event, parton)
        if (event, parton) not in split_p or (event, root) not in split_p:
            continue
        before = dphi(split_p[(event, parton)], split_p[(event, root)])
        after = dphi(momentum(r), split_p[(event, root)])
        angle_deltas.append((after - before, before, after, event, parton, r["note"]))

    out.append("### Mode-E kick angular summary")
    out.append(f"- Mode-E recursive kick records: {len(kicks)}")
    out.append(f"- coherent recursive kicks: {len(coherent_kicks)}")
    out.append(f"- resolving recursive kicks: {len(resolving_kicks)}")
    if angle_deltas:
        avg = sum(x[0] for x in angle_deltas) / len(angle_deltas)
        max_abs = max(angle_deltas, key=lambda x: abs(x[0]))
        out.append(f"- average DeltaPhi(after-before) relative to shower root: {avg:.5g}")
        out.append(
            "- largest absolute angular change: "
            f"event {max_abs[3]}, parton {max_abs[4]}, {max_abs[5]}, "
            f"before={max_abs[1]:.5g}, after={max_abs[2]:.5g}, delta={max_abs[0]:.5g}"
        )
    out.append("")

    closures = [r for r in rows if r["record_type"] == "recursive_opening_closure"]
    out.append("### Live parent-opening closure")
    out.append(f"- closure records: {len(closures)}")
    if closures:
        rel = [as_float(r["qperp"]) for r in closures]
        spatial = [extract_note_float(r["note"], "spatial_abs") or 0.0 for r in closures]
        energy = [extract_note_float(r["note"], "energy_abs") or 0.0 for r in closures]
        out.append(f"- average relative spatial residual: {sum(rel)/len(rel):.5g}")
        out.append(f"- max relative spatial residual: {max(rel):.5g}")
        out.append(f"- average absolute spatial residual: {sum(spatial)/len(spatial):.5g}")
        out.append(f"- max absolute spatial residual: {max(spatial):.5g}")
        out.append(f"- max absolute energy residual: {max(energy):.5g}")
    out.append("")
    return "\n".join(out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("histories", nargs="+", help="NAME=history.tsv or history.tsv")
    parser.add_argument("--max-examples", type=int, default=5)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    sections = ["# Mode E validation summary", ""]
    for item in args.histories:
        if "=" in item:
            name, path_s = item.split("=", 1)
            path = Path(path_s)
        else:
            path = Path(item)
            name = path.stem
        sections.append(summarize_case(name, load(path), args.max_examples))
    text = "\n".join(sections)
    if args.output:
        args.output.write_text(text)
    else:
        print(text)


if __name__ == "__main__":
    main()
