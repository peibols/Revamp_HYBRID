#!/usr/bin/env python3
"""Build an LRES tree/timeline dump and example animation from HYBRID history TSV.

The input is produced by HYBRID with:

  dump_hybrid_evolution_history = true

This script intentionally stays outside the HYBRID runtime. It converts the
time-ordered diagnostic records into a compact JSON description of the binary
shower tree, sibling-pair resolution intervals, Moliere kicks, and medium
response records, then draws an example timeline animation.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

import matplotlib.animation as animation
import matplotlib.patches as patches
import matplotlib.pyplot as plt


@dataclass
class Node:
    parton_id: int
    parent_id: int
    daughters: list[int] = field(default_factory=list)
    formation_time: float = 0.0
    x: float = 0.0
    y: float = 0.0
    z: float = 0.0
    px: float = 0.0
    py: float = 0.0
    pz: float = 0.0
    energy: float = 0.0
    depth: int = 0
    row: float = 0.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Dump and animate the LRES tree/resolution timeline."
    )
    parser.add_argument("history_tsv", type=Path, help="HYBRID evolution history TSV")
    parser.add_argument("--event", type=int, default=0, help="event_id to visualize")
    parser.add_argument("--output-dir", type=Path, required=True, help="output directory")
    parser.add_argument("--max-time", type=float, default=None, help="optional time cut")
    parser.add_argument("--fps", type=int, default=8, help="GIF frame rate")
    parser.add_argument("--frames", type=int, default=80, help="number of animation frames")
    return parser.parse_args()


def as_int(value: str) -> int:
    return int(float(value))


def as_float(value: str) -> float:
    return float(value)


def load_history(path: Path, event_id: int) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for row in reader:
            if as_int(row["event_id"]) != event_id:
                continue
            for key in ["event_id", "parton_id", "parent_id", "d1", "d2"]:
                row[key] = as_int(row[key])
            for key in ["time", "x", "y", "z", "px", "py", "pz", "E", "qperp"]:
                row[key] = as_float(row[key])
            records.append(row)
    return records


def build_tree_dump(records: list[dict[str, Any]], event_id: int) -> dict[str, Any]:
    nodes: dict[int, Node] = {}
    opening_angles: dict[tuple[int, int], float] = {}
    kicks: list[dict[str, Any]] = []
    responses: list[dict[str, Any]] = []

    for rec in records:
        rtype = rec["record_type"]
        if rtype == "split":
            pid = rec["parton_id"]
            node = nodes.get(pid, Node(parton_id=pid, parent_id=rec["parent_id"]))
            node.parent_id = rec["parent_id"]
            node.formation_time = rec["time"]
            node.x = rec["x"]
            node.y = rec["y"]
            node.z = rec["z"]
            node.px = rec["px"]
            node.py = rec["py"]
            node.pz = rec["pz"]
            node.energy = rec["E"]
            daughters = [d for d in [rec["d1"], rec["d2"]] if d >= 0]
            if daughters:
                node.daughters = daughters
            nodes[pid] = node
        elif rtype == "opening_angle":
            d1, d2 = sorted([rec["d1"], rec["d2"]])
            if d1 >= 0 and d2 >= 0:
                opening_angles[(d1, d2)] = rec["qperp"]
        elif rtype == "moliere_kick":
            kicks.append(
                {
                    "parton_id": rec["parton_id"],
                    "parent_id": rec["parent_id"],
                    "time": rec["time"],
                    "x": rec["x"],
                    "y": rec["y"],
                    "z": rec["z"],
                    "qperp": rec["qperp"],
                    "label": rec["label"],
                    "note": rec["note"],
                }
            )
        elif rtype == "medium_response":
            responses.append(
                {
                    "time": rec["time"],
                    "x": rec["x"],
                    "y": rec["y"],
                    "z": rec["z"],
                    "px": rec["px"],
                    "py": rec["py"],
                    "pz": rec["pz"],
                    "energy": rec["E"],
                    "label": rec["label"],
                    "note": rec["note"],
                }
            )

    children = defaultdict(list)
    for node in nodes.values():
        if node.parent_id >= 0:
            children[node.parent_id].append(node.parton_id)

    for pid, child_ids in children.items():
        if pid in nodes and not nodes[pid].daughters:
            nodes[pid].daughters = sorted(child_ids)

    roots = sorted(pid for pid, node in nodes.items() if node.parent_id < 0)

    def assign_depth(pid: int, depth: int) -> None:
        if pid not in nodes:
            return
        nodes[pid].depth = max(nodes[pid].depth, depth)
        for child in sorted(set(nodes[pid].daughters + children.get(pid, []))):
            assign_depth(child, depth + 1)

    for root in roots:
        assign_depth(root, 0)

    order = sorted(nodes.values(), key=lambda n: (n.depth, n.formation_time, n.parton_id))
    for i, node in enumerate(order):
        node.row = float(len(order) - 1 - i)

    intervals = []
    for parent in nodes.values():
        daughters = sorted(set(parent.daughters + children.get(parent.parton_id, [])))
        if len(daughters) < 2:
            continue
        for i in range(len(daughters)):
            for j in range(i + 1, len(daughters)):
                d1, d2 = daughters[i], daughters[j]
                if d1 not in nodes or d2 not in nodes:
                    continue
                start = max(parent.formation_time, nodes[d1].formation_time, nodes[d2].formation_time)
                end = max(start, max(nodes[d1].formation_time, nodes[d2].formation_time))
                intervals.append(
                    {
                        "parent_id": parent.parton_id,
                        "daughters": [d1, d2],
                        "start_time": start,
                        "end_time": end,
                        "opening_angle": opening_angles.get(tuple(sorted([d1, d2])), None),
                    }
                )

    return {
        "event_id": event_id,
        "description": "LRES tree/timeline dump derived from HYBRID evolution history TSV",
        "nodes": [asdict(nodes[pid]) for pid in sorted(nodes)],
        "edges": [
            {"parent_id": node.parent_id, "child_id": node.parton_id}
            for node in nodes.values()
            if node.parent_id >= 0
        ],
        "sibling_pair_intervals": intervals,
        "moliere_kicks": kicks,
        "medium_response": responses,
    }


def write_dump(dump: dict[str, Any], path: Path) -> None:
    path.write_text(json.dumps(dump, indent=2, sort_keys=True) + "\n")


def draw_frame(ax: plt.Axes, dump: dict[str, Any], cursor_time: float, max_time: float) -> None:
    ax.clear()
    nodes = {node["parton_id"]: node for node in dump["nodes"]}
    rows = {pid: node["row"] for pid, node in nodes.items()}

    ax.set_xlim(-0.05 * max_time, max_time * 1.05)
    ax.set_ylim(-1.0, max(rows.values(), default=1.0) + 1.5)
    ax.set_xlabel(r"time $\tau$ [fm/$c$]")
    ax.set_ylabel("parton tree row")
    ax.set_title("LRES tree and resolution timeline")
    ax.grid(True, axis="x", alpha=0.2)

    for interval in dump["sibling_pair_intervals"]:
        start = interval["start_time"]
        end = min(interval["end_time"], max_time)
        if start > max_time:
            continue
        d1, d2 = interval["daughters"]
        y0 = min(rows.get(d1, 0), rows.get(d2, 0)) - 0.25
        height = abs(rows.get(d1, 0) - rows.get(d2, 0)) + 0.5
        rect = patches.Rectangle(
            (start, y0),
            max(end - start, 0.05 * max_time),
            height,
            facecolor="0.88",
            edgecolor="0.55",
            lw=1.0,
            alpha=0.45,
            zorder=0,
        )
        ax.add_patch(rect)
        ax.text(start + 0.02 * max_time, y0 + height - 0.18, r"$\Delta R<L_{\rm res}$", fontsize=8, color="0.25")

    for edge in dump["edges"]:
        parent = nodes.get(edge["parent_id"])
        child = nodes.get(edge["child_id"])
        if parent is None or child is None:
            continue
        split_time = child["formation_time"]
        if split_time > max_time:
            continue
        ax.plot([parent["formation_time"], split_time], [parent["row"], parent["row"]], color="black", lw=1.8)
        ax.plot([split_time, split_time], [parent["row"], child["row"]], color="black", lw=1.0, alpha=0.7)

    for node in nodes.values():
        t0 = node["formation_time"]
        if t0 > max_time:
            continue
        color = "tab:blue" if node["parent_id"] >= 0 else "black"
        ax.plot([t0, max_time], [node["row"], node["row"]], color=color, lw=1.2, alpha=0.35)
        ax.scatter([t0], [node["row"]], s=max(20.0, min(120.0, node["energy"] / 8.0)), color=color, zorder=4)
        ax.text(t0, node["row"] + 0.16, f"id {node['parton_id']}, E={node['energy']:.1f}", fontsize=7)

    for kick in dump["moliere_kicks"]:
        t = kick["time"]
        if t > max_time:
            continue
        pid = kick["parton_id"]
        y = rows.get(pid, max(rows.values(), default=0) + 0.5)
        length = 0.15 + min(kick["qperp"], 50.0) / 70.0
        ax.annotate(
            "",
            xy=(min(t + 0.06 * max_time, max_time), y + length),
            xytext=(t, y),
            arrowprops={"arrowstyle": "->", "color": "red", "lw": 1.4, "alpha": 0.85},
            zorder=5,
        )
        ax.text(t, y + length + 0.08, rf"$q_\perp={kick['qperp']:.1f}$", fontsize=7, color="red")

    for response in dump["medium_response"]:
        t = response["time"]
        if t > max_time:
            continue
        ax.scatter([t], [-0.45], s=45, marker="D", color="tab:green", alpha=0.75, zorder=5)
        ax.text(t, -0.75, response["note"], fontsize=7, color="tab:green", rotation=30)

    ax.axvline(cursor_time, color="tab:orange", lw=2.0, alpha=0.8)
    ax.text(cursor_time, ax.get_ylim()[1] - 0.35, rf"$\tau={cursor_time:.2f}$", color="tab:orange", fontsize=10)

    legend_handles = [
        patches.Patch(facecolor="0.88", edgecolor="0.55", alpha=0.45, label=r"unresolved by $L_{\rm res}$"),
        plt.Line2D([0], [0], color="black", lw=1.8, label="coherent/tree branch"),
        plt.Line2D([0], [0], color="tab:blue", lw=1.8, label="resolved daughter timeline"),
        plt.Line2D([0], [0], color="red", lw=1.8, label="Moliere kick"),
        plt.Line2D([0], [0], marker="D", color="tab:green", lw=0, label="medium response"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=8, frameon=True)


def make_animation(dump: dict[str, Any], output_dir: Path, max_time: float, fps: int, frames: int) -> None:
    fig, ax = plt.subplots(figsize=(11, 6.2), constrained_layout=True)
    times = [max_time * i / max(frames - 1, 1) for i in range(frames)]

    def update(i: int) -> None:
        draw_frame(ax, dump, times[i], max_time)

    anim = animation.FuncAnimation(fig, update, frames=len(times), interval=1000 / fps)
    gif_path = output_dir / f"lres_tree_timeline_event{dump['event_id']}.gif"
    png_path = output_dir / f"lres_tree_timeline_event{dump['event_id']}.png"
    update(len(times) - 1)
    fig.savefig(png_path, dpi=160)
    anim.save(gif_path, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records = load_history(args.history_tsv, args.event)
    if not records:
        raise RuntimeError(f"No records found for event_id={args.event} in {args.history_tsv}")

    dump = build_tree_dump(records, args.event)
    dump_path = args.output_dir / f"lres_tree_timeline_event{args.event}.json"
    write_dump(dump, dump_path)

    all_times = [node["formation_time"] for node in dump["nodes"]]
    all_times += [kick["time"] for kick in dump["moliere_kicks"]]
    all_times += [response["time"] for response in dump["medium_response"]]
    max_time = args.max_time if args.max_time is not None else max(all_times + [1.0])
    make_animation(dump, args.output_dir, max_time, args.fps, args.frames)
    print(f"Wrote {dump_path}")
    print(f"Wrote {args.output_dir / f'lres_tree_timeline_event{args.event}.png'}")
    print(f"Wrote {args.output_dir / f'lres_tree_timeline_event{args.event}.gif'}")


if __name__ == "__main__":
    main()
