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
from collections import Counter, defaultdict
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

    kicks = [r for r in rows if r["record_type"] == "moliere_kick" and r["note"].startswith("modeE_")]
    coherent_kicks = [
        r
        for r in kicks
        if "coherent" in r["note"] or "active_coherent_source" in r["note"]
    ]
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

    failed_daughters = [r for r in rows if r["record_type"] == "recursive_failed_daughter_veto"]
    parent_tests = [r for r in rows if r["record_type"] == "recursive_coherent_resample_test"]
    parent_accepts = [r for r in parent_tests if "accept_coherent" in r["note"]]
    parent_vetoes = [r for r in parent_tests if "veto_resolving" in r["note"]]
    parent_exhausted = [
        r for r in rows if r["record_type"] == "recursive_coherent_resample_exhausted"
    ]
    neutral_skips = [r for r in rows if r["record_type"] == "recursive_coherent_resample_skip"]
    resampled_kicks = [r for r in kicks if "resampled_coherent_parent" in r["note"]]
    direct_coherent_kicks = [r for r in kicks if "active_coherent_source" in r["note"]]

    accepted_keys = {
        (as_int(r["event_id"]), as_int(r["parton_id"]), r["time"])
        for r in parent_accepts
    }
    resampled_keys = {
        (as_int(r["event_id"]), as_int(r["parton_id"]), r["time"])
        for r in resampled_kicks
    }
    unmatched_accepts = accepted_keys - resampled_keys
    unmatched_kicks = resampled_keys - accepted_keys

    chosen_frontier_candidates = (
        len(failed_daughters) + len(resolving_kicks) + len(direct_coherent_kicks)
    )
    classified_parent_candidates = len(parent_accepts) + len(parent_vetoes)
    classified_all_candidates = chosen_frontier_candidates + classified_parent_candidates

    out.append("### Dani failed-probe veto/resampling contract")
    out.append(f"- failed daughter candidates vetoed: {len(failed_daughters)}")
    out.append(f"- coherent-source candidates inspected: {len(parent_tests)}")
    out.append(f"- coherent-source candidates accepted as unresolved: {len(parent_accepts)}")
    out.append(f"- coherent-source candidates vetoed as resolving: {len(parent_vetoes)}")
    out.append(f"- coherent searches reaching the LRES boundary without acceptance: {len(parent_exhausted)}")
    out.append(f"- color-neutral coherent sources skipped: {len(neutral_skips)}")
    out.append(f"- accepted coherent-source kicks committed: {len(resampled_kicks)}")
    out.append(
        "- parent candidate classification closure: "
        f"{classified_parent_candidates}/{len(parent_tests)} "
        f"({'PASS' if classified_parent_candidates == len(parent_tests) else 'FAIL'})"
    )
    out.append(
        "- accepted-parent-test to committed-kick matching: "
        f"{'PASS' if not unmatched_accepts and not unmatched_kicks else 'FAIL'}"
    )
    out.append(
        "- history-classified hard candidates "
        "(chosen frontier + coherent resamples): "
        f"{classified_all_candidates}"
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


def validate_dani_contract(
    rows: List[Dict[str, str]],
    require_failed_veto: bool,
    require_parent_accept: bool,
    require_parent_veto: bool,
    require_neutral_skip: bool,
) -> List[str]:
    errors: List[str] = []
    failed = [r for r in rows if r["record_type"] == "recursive_failed_daughter_veto"]
    tests = [r for r in rows if r["record_type"] == "recursive_coherent_resample_test"]
    accepts = [r for r in tests if "accept_coherent" in r["note"]]
    vetoes = [r for r in tests if "veto_resolving" in r["note"]]
    exhausted = [r for r in rows if r["record_type"] == "recursive_coherent_resample_exhausted"]
    neutral = [r for r in rows if r["record_type"] == "recursive_coherent_resample_skip"]
    kicks = [
        r
        for r in rows
        if r["record_type"] == "moliere_kick"
        and "resampled_coherent_parent" in r["note"]
    ]
    old_mapped = [
        r
        for r in rows
        if r["record_type"] == "moliere_kick"
        and r["note"] == "modeE_recursive_coherent_kick"
    ]

    if len(tests) != len(accepts) + len(vetoes):
        errors.append("not every coherent-source candidate is classified as accept or veto")
    if len(failed) != len(accepts) + len(exhausted) + len(neutral):
        errors.append("failed daughter probes do not close against accept/exhaust/neutral outcomes")

    def key(row: Dict[str, str]) -> Tuple[int, int, str]:
        return (as_int(row["event_id"]), as_int(row["parton_id"]), row["time"])

    if Counter(map(key, accepts)) != Counter(map(key, kicks)):
        errors.append("accepted coherent-source candidates do not match committed parent kicks")
    if old_mapped:
        errors.append("obsolete daughter-sampled coherent kick records remain")
    if require_failed_veto and not failed:
        errors.append("required failed-daughter veto path was not exercised")
    if require_parent_accept and not accepts:
        errors.append("required coherent-source acceptance path was not exercised")
    if require_parent_veto and not vetoes:
        errors.append("required resolving parent-candidate veto path was not exercised")
    if require_neutral_skip and not neutral:
        errors.append("required color-neutral coherent-source skip was not exercised")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("histories", nargs="+", help="NAME=history.tsv or history.tsv")
    parser.add_argument("--max-examples", type=int, default=5)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--strict", action="store_true", help="exit nonzero on contract failure")
    parser.add_argument("--require-failed-veto", action="store_true")
    parser.add_argument("--require-parent-accept", action="store_true")
    parser.add_argument("--require-parent-veto", action="store_true")
    parser.add_argument("--require-neutral-skip", action="store_true")
    args = parser.parse_args()

    sections = ["# Mode E validation summary", ""]
    validation_errors: List[str] = []
    for item in args.histories:
        if "=" in item:
            name, path_s = item.split("=", 1)
            path = Path(path_s)
        else:
            path = Path(item)
            name = path.stem
        rows = load(path)
        sections.append(summarize_case(name, rows, args.max_examples))
        if args.strict:
            validation_errors.extend(
                f"{name}: {error}"
                for error in validate_dani_contract(
                    rows,
                    args.require_failed_veto,
                    args.require_parent_accept,
                    args.require_parent_veto,
                    args.require_neutral_skip,
                )
            )
    text = "\n".join(sections)
    if args.output:
        args.output.write_text(text)
    else:
        print(text)
    if validation_errors:
        for error in validation_errors:
            print(f"ERROR: {error}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
