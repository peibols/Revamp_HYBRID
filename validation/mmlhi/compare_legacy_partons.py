#!/usr/bin/env python3
"""Compare same-tree legacy-heavy and MMLHI parton outputs."""

from __future__ import annotations

import argparse
import math
import statistics
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Parton:
    px: float
    py: float
    pz: float
    mass: float
    pdg_id: int
    status: int

    @property
    def pt(self) -> float:
        return math.hypot(self.px, self.py)

    @property
    def momentum(self) -> float:
        return math.sqrt(self.px * self.px + self.py * self.py + self.pz * self.pz)

    @property
    def energy(self) -> float:
        return math.sqrt(self.momentum * self.momentum + self.mass * self.mass)

    @property
    def eta(self) -> float:
        if self.pt == 0.:
            return math.copysign(math.inf, self.pz)
        return math.asinh(self.pz / self.pt)


@dataclass(frozen=True)
class Event:
    event_id: int
    metadata: tuple[float, ...]
    partons: tuple[Parton, ...]


def parse(path: Path) -> list[Event]:
    events: list[Event] = []
    event_id: int | None = None
    metadata: tuple[float, ...] = ()
    partons: list[Parton] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        fields = raw_line.split()
        if not fields:
            continue
        if fields[:2] == ["#", "event"]:
            if event_id is not None:
                raise RuntimeError(f"{path}: event {event_id} has no end marker")
            event_id = int(fields[2])
            metadata = ()
            partons = []
        elif fields[0] == "weight":
            metadata = tuple(float(fields[index]) for index in (1, 3, 5, 7))
        elif fields[0] == "end":
            if event_id is None:
                raise RuntimeError(f"{path}: end marker outside an event")
            events.append(Event(event_id, metadata, tuple(partons)))
            event_id = None
        elif len(fields) == 6:
            partons.append(
                Parton(
                    px=float(fields[0]),
                    py=float(fields[1]),
                    pz=float(fields[2]),
                    mass=float(fields[3]),
                    pdg_id=int(fields[4]),
                    status=int(fields[5]),
                )
            )
    if event_id is not None:
        raise RuntimeError(f"{path}: final event has no end marker")
    return events


def relative_difference(reference: float, candidate: float) -> float:
    return (candidate - reference) / reference if reference != 0. else 0.


def opening_angle(first: Parton, second: Parton) -> float:
    denominator = first.momentum * second.momentum
    if denominator == 0.:
        return 0.
    cosine = (first.px * second.px + first.py * second.py + first.pz * second.pz) / denominator
    return math.acos(max(-1., min(1., cosine)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("legacy", type=Path)
    parser.add_argument("mmlhi", type=Path)
    parser.add_argument("--max-abs-eta", type=float)
    args = parser.parse_args()

    legacy_events = parse(args.legacy)
    mmlhi_events = parse(args.mmlhi)
    if len(legacy_events) != len(mmlhi_events):
        raise RuntimeError("event counts differ")

    initial_max_abs_difference = 0.
    pt_relative_differences: list[float] = []
    pt_absolute_differences: list[float] = []
    energy_relative_differences: list[float] = []
    angular_differences: list[float] = []
    legacy_pts: list[float] = []
    mmlhi_pts: list[float] = []
    comparison_labels: list[tuple[int, int]] = []

    for legacy_event, mmlhi_event in zip(legacy_events, mmlhi_events):
        if legacy_event.event_id != mmlhi_event.event_id:
            raise RuntimeError("event identifiers differ")
        if legacy_event.metadata != mmlhi_event.metadata:
            raise RuntimeError(f"event {legacy_event.event_id}: metadata differ")

        legacy_initial = [parton for parton in legacy_event.partons if parton.status == -2]
        mmlhi_initial = [parton for parton in mmlhi_event.partons if parton.status == -2]
        if len(legacy_initial) != len(mmlhi_initial):
            raise RuntimeError(f"event {legacy_event.event_id}: initial hard-parton counts differ")
        for legacy_parton, mmlhi_parton in zip(legacy_initial, mmlhi_initial):
            if (legacy_parton.pdg_id, legacy_parton.status) != (
                mmlhi_parton.pdg_id,
                mmlhi_parton.status,
            ):
                raise RuntimeError(f"event {legacy_event.event_id}: initial identities differ")
            initial_max_abs_difference = max(
                initial_max_abs_difference,
                *(abs(left - right) for left, right in zip(
                    (legacy_parton.px, legacy_parton.py, legacy_parton.pz, legacy_parton.mass),
                    (mmlhi_parton.px, mmlhi_parton.py, mmlhi_parton.pz, mmlhi_parton.mass),
                )),
            )

        legacy_heavy = [
            parton for parton in legacy_event.partons
            if parton.status == 0 and abs(parton.pdg_id) in (4, 5)
        ]
        mmlhi_heavy = [
            parton for parton in mmlhi_event.partons
            if parton.status == 0 and abs(parton.pdg_id) in (4, 5)
        ]
        if [parton.pdg_id for parton in legacy_heavy] != [parton.pdg_id for parton in mmlhi_heavy]:
            raise RuntimeError(f"event {legacy_event.event_id}: final heavy identities differ")
        for legacy_parton, mmlhi_parton in zip(legacy_heavy, mmlhi_heavy):
            if args.max_abs_eta is not None and (
                abs(legacy_parton.eta) >= args.max_abs_eta or
                abs(mmlhi_parton.eta) >= args.max_abs_eta
            ):
                continue
            comparison_labels.append((legacy_event.event_id, legacy_parton.pdg_id))
            legacy_pts.append(legacy_parton.pt)
            mmlhi_pts.append(mmlhi_parton.pt)
            pt_relative_differences.append(relative_difference(legacy_parton.pt, mmlhi_parton.pt))
            pt_absolute_differences.append(mmlhi_parton.pt - legacy_parton.pt)
            energy_relative_differences.append(
                relative_difference(legacy_parton.energy, mmlhi_parton.energy)
            )
            angular_differences.append(opening_angle(legacy_parton, mmlhi_parton))

    def rms(values: list[float]) -> float:
        return math.sqrt(statistics.fmean(value * value for value in values)) if values else 0.

    def percentile(values: list[float], fraction: float) -> float:
        ordered = sorted(values)
        position = fraction * (len(ordered) - 1)
        lower = math.floor(position)
        upper = math.ceil(position)
        if lower == upper:
            return ordered[lower]
        return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])

    if not pt_relative_differences:
        raise RuntimeError("no final heavy partons pass the requested acceptance")
    absolute_fractional_pt_differences = list(map(abs, pt_relative_differences))
    maximum_index = max(
        range(len(absolute_fractional_pt_differences)),
        key=absolute_fractional_pt_differences.__getitem__,
    )
    maximum_event, maximum_pdg_id = comparison_labels[maximum_index]

    metrics: tuple[tuple[str, int | float], ...] = (
        ("events", len(legacy_events)),
        ("max_abs_eta", args.max_abs_eta if args.max_abs_eta is not None else math.inf),
        ("matched_final_heavy", len(pt_relative_differences)),
        ("initial_max_abs_difference", initial_max_abs_difference),
        ("legacy_mean_pt", statistics.fmean(legacy_pts)),
        ("mmlhi_mean_pt", statistics.fmean(mmlhi_pts)),
        ("mean_fractional_pt_difference", statistics.fmean(pt_relative_differences)),
        ("rms_fractional_pt_difference", rms(pt_relative_differences)),
        ("median_abs_fractional_pt_difference", statistics.median(absolute_fractional_pt_differences)),
        ("p90_abs_fractional_pt_difference", percentile(absolute_fractional_pt_differences, 0.9)),
        ("max_abs_fractional_pt_difference", max(absolute_fractional_pt_differences)),
        ("max_difference_event", maximum_event),
        ("max_difference_pdg_id", maximum_pdg_id),
        ("max_difference_legacy_pt_gev", legacy_pts[maximum_index]),
        ("max_difference_mmlhi_pt_gev", mmlhi_pts[maximum_index]),
        ("mean_abs_pt_difference_gev", statistics.fmean(map(abs, pt_absolute_differences))),
        ("max_abs_pt_difference_gev", max(map(abs, pt_absolute_differences))),
        ("mean_fractional_energy_difference", statistics.fmean(energy_relative_differences)),
        ("rms_fractional_energy_difference", rms(energy_relative_differences)),
        ("mean_direction_difference_rad", statistics.fmean(angular_differences)),
        ("max_direction_difference_rad", max(angular_differences)),
    )
    print("metric\tvalue")
    for name, value in metrics:
        print(f"{name}\t{value:.12g}" if isinstance(value, float) else f"{name}\t{value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
