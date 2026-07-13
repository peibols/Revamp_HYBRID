#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import struct
import sys
import time
from typing import Callable

CENTRALITIES = ["C0-5", "C5-10", "C10-20", "C20-30", "C30-40", "C40-50", "C50-60", "C60-70", "C70-80", "C80-90", "C90-100"]
DEFAULT_PREHYDRO_TAU_GRID = ",".join(f"{i / 100:.2f}" for i in range(1, 40)) + ",0.399"
HBARC_GEV_FM = 0.1973269804
PUBLISHED_PREHYDRO_ETA_OVER_S = 0.12
QCD_KINETIC_ATTRACTOR_C_INF = 0.87351
QCD_THREE_FLAVOR_EFFECTIVE_DEGREES = 47.5
QCD_CONFORMAL_EOS_FACTOR = math.pi**2 * QCD_THREE_FLAVOR_EFFECTIVE_DEGREES / 30.0
PUBLISHED_QCD_ATTRACTOR_SHA256 = "1bea7289d3dc8ed95819eaa86cf4c489442a054c14aae47eff010cf45155eba0"
PUBLISHED_QCD_ATTRACTOR_DOI = "10.4119/unibi/2939684"
PLAN_B_REFERENCE = "arXiv:2509.19430v2 Eqs. (2)-(4)"
DEFAULT_NO_PREHYDRO_ALPHA = 0.37
DEFAULT_PREHYDRO_ALPHA = 0.37
DEFAULT_BROADENING_K = 15.0


def b(v: bool) -> str:
    return "true" if v else "false"


NPDF_KEYS = {
    "PDF:useHardNPDFB",
    "PDF:nPDFSetB",
    "PDF:nPDFBeamB",
    "PDF:useHardNPDFA",
    "PDF:nPDFSetA",
    "PDF:nPDFBeamA",
}


def write_pythia_card(template: Path, out: Path, pthat_min: float, pthat_max: float, pdf_mode: str, lhapdf_set: str) -> None:
    lines = []
    replacements = {
        "Beams:eCM": "Beams:eCM = 5360.",
        "PhaseSpace:bias2Selection": "PhaseSpace:bias2Selection = on",
        "PhaseSpace:bias2SelectionPow": "PhaseSpace:bias2SelectionPow = 4.",
        "PhaseSpace:bias2SelectionRef": "PhaseSpace:bias2SelectionRef = 10.",
        "PhaseSpace:pTHatMin": f"PhaseSpace:pTHatMin = {pthat_min:g}",
        "PhaseSpace:pTHatMax": f"PhaseSpace:pTHatMax = {pthat_max:g}",
    }
    if pdf_mode == "native_npdf":
        replacements.update({
            "PDF:pSet": "PDF:pSet = 13",
            "PDF:useHardNPDFA": "PDF:useHardNPDFA = on",
            "PDF:useHardNPDFB": "PDF:useHardNPDFB = on",
            "PDF:nPDFSetA": "PDF:nPDFSetA = 1",
            "PDF:nPDFSetB": "PDF:nPDFSetB = 1",
            "PDF:nPDFBeamA": "PDF:nPDFBeamA = 1000080160",
            "PDF:nPDFBeamB": "PDF:nPDFBeamB = 1000080160",
        })
    elif pdf_mode == "lhapdf":
        replacements["PDF:pSet"] = f"PDF:pSet = LHAPDF6:{lhapdf_set}"
    elif pdf_mode == "off":
        replacements["PDF:pSet"] = "PDF:pSet = 13"
    else:
        raise ValueError(f"Unknown PDF mode: {pdf_mode}")

    seen = set()
    for raw in template.read_text().splitlines():
        stripped = raw.strip()
        uncommented = stripped[1:].strip() if stripped.startswith("!") else stripped
        key = uncommented.split("=", 1)[0].strip() if "=" in uncommented else None
        if pdf_mode in {"off", "lhapdf"} and key in NPDF_KEYS:
            lines.append(raw if stripped.startswith("!") else f"!{raw}")
            seen.add(key)
            continue
        if key in replacements:
            lines.append(replacements[key])
            seen.add(key)
        else:
            lines.append(raw)
    for key, value in replacements.items():
        if key not in seen:
            lines.append(value)
    out.write_text("\n".join(lines) + "\n")


def write_input(
    path: Path,
    *,
    seed: int,
    events: int,
    aa: bool,
    energy_loss_alpha: float = DEFAULT_NO_PREHYDRO_ALPHA,
    broadening_k: float = DEFAULT_BROADENING_K,
    use_prehydro: bool = False,
    prehydro_file: str = "prehydro_table.tsv",
) -> None:
    lines = [
        f"njob = {seed}",
        f"seed_base = {seed}",
        f"Nev = {events}",
        "cent = 0-100",
        f"alpha = {energy_loss_alpha}",
        f"kappa = {broadening_k}",
        "tmethod = 0",
        "mode = 0",
        f"do_quench = {b(aa)}",
        f"do_wake = {b(aa)}",
        "do_source = false",
        "do_elastic = false",
        "do_lres = false",
        "do_Moliere_on_unresolved_partons = false",
        "do_Moliere_dynamic_unresolved_resolution = false",
        "do_Moliere_dynamic_daughter_unresolved_resolution = false",
        f"ebe_hydro = {1 if aa else 0}",
        "hadro_type = 0",
        "compat_moliere_legacy_hydro = false",
        f"use_prehydro = {b(use_prehydro)}",
        f"prehydro_file = {prehydro_file}",
        "output_base = HYBRID",
    ]
    path.write_text("\n".join(lines) + "\n")


def run_main(run_dir: Path, main_bin: Path, timeout: int) -> int:
    with (run_dir / "stdout.log").open("w") as out:
        proc = subprocess.run([str(main_bin), "hybrid_input.dat"], cwd=run_dir, stdout=out, stderr=subprocess.STDOUT, timeout=timeout)
    return proc.returncode


def parse_tau_grid(raw: str) -> list[float]:
    values = [float(item.strip()) for item in raw.split(",") if item.strip()]
    if not values:
        raise ValueError("prehydro tau grid is empty")
    if values != sorted(values):
        raise ValueError(f"prehydro tau grid must be sorted: {values}")
    return values


def parse_hydro_header(path: Path) -> dict[str, float | int]:
    with path.open("rb") as stream:
        raw = stream.read(64)
    if len(raw) != 64:
        raise RuntimeError(f"{path} is too small to contain a 16-float hydro header")
    header = struct.unpack("<16f", raw)
    return {
        "tau0_fm": float(header[0]),
        "dtau_fm": float(header[1]),
        "nx": int(header[2]),
        "dx_fm": float(header[3]),
        "x_extent_fm": abs(float(header[4])),
        "ny": int(header[5]),
        "dy_fm": float(header[6]),
        "y_extent_fm": abs(float(header[7])),
        "n_fields": int(header[15]),
    }


def first_hydro_slice(path: Path) -> tuple[dict[str, float | int], list[dict[str, float | int]]]:
    header = parse_hydro_header(path)
    n_fields = int(header["n_fields"])
    if n_fields < 11:
        raise RuntimeError(f"{path}: expected at least 11 fields per hydro cell, got {n_fields}")
    x0 = -float(header["x_extent_fm"])
    y0 = -float(header["y_extent_fm"])
    dx = float(header["dx_fm"])
    dy = float(header["dy_fm"])
    rows: list[dict[str, float | int]] = []
    record_size = 4 * n_fields
    with path.open("rb") as stream:
        stream.seek(64)
        while True:
            raw = stream.read(record_size)
            if not raw:
                break
            if len(raw) != record_size:
                raise RuntimeError(f"{path}: truncated hydro record")
            fields = struct.unpack("<" + "f" * n_fields, raw)
            it = int(fields[0])
            if it != 0:
                continue
            ix = int(fields[1])
            iy = int(fields[2])
            temperature = float(fields[6])
            ux = float(fields[8])
            uy = float(fields[9])
            uz = float(fields[10])
            gamma = math.sqrt(1.0 + ux * ux + uy * uy + uz * uz)
            rows.append(
                {
                    "ix": ix,
                    "iy": iy,
                    "x_fm": x0 + ix * dx,
                    "y_fm": y0 + iy * dy,
                    "temperature_hyd": temperature,
                    "vx_hyd": ux / gamma,
                    "vy_hyd": uy / gamma,
                }
            )
    return header, rows


def load_attractor_table(path: Path) -> tuple[str, Callable[[float], float]]:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != PUBLISHED_QCD_ATTRACTOR_SHA256:
        raise ValueError(
            f"{path}: expected published attractor SHA256 "
            f"{PUBLISHED_QCD_ATTRACTOR_SHA256}, got {digest}"
        )
    points: list[tuple[float, float]] = []
    with path.open("r", encoding="utf-8", errors="replace") as stream:
        for raw in stream:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = re.split(r"[\s,]+", line)
            if len(fields) < 2:
                continue
            try:
                omega = float(fields[0])
                value = float(fields[1])
            except ValueError:
                continue
            if not math.isfinite(omega) or not math.isfinite(value):
                raise ValueError(f"{path}: non-finite attractor row: {line}")
            if omega <= 0.0 or value <= 0.0:
                raise ValueError(f"{path}: omega and E(omega) must be positive: {line}")
            points.append((omega, value))
    if len(points) < 2:
        raise ValueError(f"{path} must contain at least two numeric omega,E rows")
    points = sorted(points)
    for previous, current in zip(points, points[1:]):
        if current[0] <= previous[0]:
            raise ValueError(f"{path}: omega values must be strictly increasing")

    def interpolate(omega: float) -> float:
        if omega < points[0][0] or omega > points[-1][0]:
            raise ValueError(
                f"omega={omega:.8g} lies outside the published attractor range "
                f"[{points[0][0]:.8g}, {points[-1][0]:.8g}]"
            )
        if omega == points[0][0]:
            return points[0][1]
        if omega == points[-1][0]:
            return points[-1][1]
        lo = 0
        hi = len(points) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if points[mid][0] <= omega:
                lo = mid
            else:
                hi = mid
        x0, y0 = points[lo]
        x1, y1 = points[hi]
        frac = (omega - x0) / (x1 - x0)
        return y0 + frac * (y1 - y0)

    return (
        f"published_qcd_lambda10_Cinf0p87;sha256={digest};file={path}",
        interpolate,
    )


def attractor_temperature(
    *,
    tau: float,
    tau_hyd: float,
    temperature_hyd: float,
    eta_over_s: float,
    attractor: Callable[[float], float],
    viscous_anchor: bool,
) -> tuple[float, float, float]:
    if tau <= 0:
        raise ValueError(f"tau must be positive, got {tau}")
    if temperature_hyd <= 0:
        return 0.0, 0.0, 0.0

    tau_natural = tau / HBARC_GEV_FM
    tau_hyd_natural = tau_hyd / HBARC_GEV_FM
    anchor_temperature = temperature_hyd
    if viscous_anchor:
        anchor_temperature += 2.0 * eta_over_s / (3.0 * tau_hyd_natural)
    anchor = (tau_hyd_natural ** (1.0 / 3.0)) * anchor_temperature
    temperature = anchor / (tau_natural ** (1.0 / 3.0))
    omega = 0.0
    e_value = 1.0
    if eta_over_s <= 0:
        return temperature, omega, e_value

    converged = False
    for _ in range(100):
        omega = tau_natural * temperature / (4.0 * math.pi * eta_over_s)
        e_value = max(0.0, attractor(omega))
        updated = anchor * (e_value ** 0.25) / (tau_natural ** (1.0 / 3.0))
        if abs(updated - temperature) < 1e-10 * max(1.0, temperature):
            temperature = updated
            converged = True
            break
        temperature = 0.5 * temperature + 0.5 * updated
    if not converged:
        raise RuntimeError(
            f"attractor temperature did not converge at tau={tau}, "
            f"T_hyd={temperature_hyd}"
        )
    omega = tau_natural * temperature / (4.0 * math.pi * eta_over_s)
    e_value = max(0.0, attractor(omega))
    residual = temperature - anchor * (e_value ** 0.25) / (tau_natural ** (1.0 / 3.0))
    if abs(residual) > 1e-9 * max(1.0, temperature):
        raise RuntimeError(
            f"attractor temperature residual is too large at tau={tau}: {residual}"
        )
    return temperature, omega, e_value


def event_id_from_readme(hydro_dir: Path, fallback: int) -> int:
    readme = hydro_dir / "README_staged_event.txt"
    if not readme.exists():
        return fallback
    match = re.search(r"^event_id\s*=\s*(\d+)\s*$", readme.read_text(errors="replace"), re.MULTILINE)
    return int(match.group(1)) if match else fallback


AA_TASK_FIELDS = {
    "task_id",
    "hard_seed",
    "milestone_block",
    "hydro_slot",
    "hydro_event_id",
    "hydro_ncoll",
    "hydro_dir",
    "hydro_payload_key",
    "hydro_payload_sha256",
}


def load_aa_task_assignment(path: Path, task_id: int) -> dict[str, str | int]:
    if not path.is_file():
        raise FileNotFoundError(path)
    selected: dict[str, str] | None = None
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = AA_TASK_FIELDS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: task manifest is missing {sorted(missing)}")
        for row in reader:
            try:
                row_task_id = int(row["task_id"])
            except (TypeError, ValueError) as error:
                raise ValueError(f"{path}: invalid task_id row") from error
            if row_task_id == task_id:
                if selected is not None:
                    raise ValueError(f"{path}: duplicate task_id {task_id}")
                selected = row
    if selected is None:
        raise ValueError(f"{path}: no assignment for task_id {task_id}")

    integer_fields = (
        "task_id",
        "hard_seed",
        "milestone_block",
        "hydro_slot",
        "hydro_event_id",
        "hydro_ncoll",
    )
    assignment: dict[str, str | int] = dict(selected)
    for field in integer_fields:
        try:
            assignment[field] = int(selected[field])
        except (TypeError, ValueError) as error:
            raise ValueError(f"{path}: invalid integer {field} for task {task_id}") from error
    hydro_dir = str(assignment["hydro_dir"])
    expected_dir = f"C0-5_event_{int(assignment['hydro_event_id']):05d}"
    if hydro_dir != expected_dir or Path(hydro_dir).name != hydro_dir:
        raise ValueError(
            f"{path}: task {task_id} hydro_dir={hydro_dir!r}, expected {expected_dir!r}"
        )
    payload_sha = str(assignment["hydro_payload_sha256"])
    if not re.fullmatch(r"[0-9a-f]{64}", payload_sha):
        raise ValueError(f"{path}: task {task_id} has invalid hydro payload SHA256")
    if int(assignment["hydro_ncoll"]) <= 0:
        raise ValueError(f"{path}: task {task_id} has nonpositive Ncoll")
    return assignment


def read_staged_hydro_metadata(hydro_dir: Path) -> dict[str, int]:
    readme = hydro_dir / "README_staged_event.txt"
    if not readme.is_file():
        raise FileNotFoundError(readme)
    metadata: dict[str, str] = {}
    for raw in readme.read_text(errors="strict").splitlines():
        if "=" not in raw or raw.lstrip().startswith("#"):
            continue
        key, value = (part.strip() for part in raw.split("=", 1))
        metadata[key] = value
    required = ("hydro_slot", "event_id", "ncoll_positions")
    missing = [key for key in required if key not in metadata]
    if missing:
        raise ValueError(f"{readme}: missing {missing}")
    parsed = {key: int(metadata[key]) for key in required}
    ncoll_path = hydro_dir / "NcollList.dat"
    if not ncoll_path.is_file():
        raise FileNotFoundError(ncoll_path)
    actual_ncoll = sum(
        1
        for line in ncoll_path.read_text(errors="strict").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    )
    if actual_ncoll != parsed["ncoll_positions"]:
        raise ValueError(
            f"{ncoll_path}: contains {actual_ncoll} positions, "
            f"README records {parsed['ncoll_positions']}"
        )
    return parsed


def alpha_tag(value: float) -> str:
    return f"{value:.12g}".replace("-", "m").replace(".", "p").replace("+", "")


def build_aa_variants(
    *,
    base_run_dir: Path,
    task_id: int,
    run_prehydro_pair: bool,
    no_prehydro_alpha: float,
    prehydro_alpha: float,
    additional_prehydro_alphas: list[float],
    run_prehydro_only: bool = False,
) -> list[tuple[str, Path, bool, float]]:
    if run_prehydro_pair and run_prehydro_only:
        raise ValueError("prehydro pair and prehydro-only modes are exclusive")
    if run_prehydro_only:
        if additional_prehydro_alphas:
            raise ValueError(
                "additional prehydro alphas require --run-prehydro-pair"
            )
        return [
            (
                "with_prehydro",
                base_run_dir.with_name(f"task_{task_id:05d}_prehydro"),
                True,
                prehydro_alpha,
            )
        ]

    variants: list[tuple[str, Path, bool, float]] = [
        ("no_prehydro", base_run_dir, False, no_prehydro_alpha)
    ]
    if not run_prehydro_pair:
        if additional_prehydro_alphas:
            raise ValueError(
                "additional prehydro alphas require --run-prehydro-pair"
            )
        return variants

    variants.append(
        (
            "with_prehydro",
            base_run_dir.with_name(f"task_{task_id:05d}_prehydro"),
            True,
            prehydro_alpha,
        )
    )
    seen = {prehydro_alpha}
    for alpha in additional_prehydro_alphas:
        if alpha in seen:
            raise ValueError(f"duplicate prehydro alpha {alpha}")
        seen.add(alpha)
        tag = alpha_tag(alpha)
        variants.append(
            (
                f"with_prehydro_alpha_{tag}",
                base_run_dir.with_name(
                    f"task_{task_id:05d}_prehydro_alpha_{tag}"
                ),
                True,
                alpha,
            )
        )
    return variants


def generate_planb_prehydro(
    *,
    reference_hydro: Path,
    output_path: Path,
    event_id: int,
    tau_grid: list[float],
    tau_min: float,
    eta_over_s: float,
    eos_factor: float,
    attractor_label: str,
    attractor: Callable[[float], float],
    viscous_anchor: bool,
) -> None:
    header, cells = first_hydro_slice(reference_hydro)
    tau_hyd = float(header["tau0_fm"])
    columns = [
        "event_id",
        "tau_fm",
        "ix",
        "iy",
        "x_fm",
        "y_fm",
        "temperature_eff",
        "energy_density_eff",
        "vx",
        "vy",
        "valid",
        "temperature_hyd_tau0",
        "vx_hyd_tau0",
        "vy_hyd_tau0",
        "attractor_omega",
        "attractor_E",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as stream:
        stream.write("# arXiv:2509.19430v2 hydro-anchored pre-equilibrium table\n")
        stream.write(f"# reference_hydro = {reference_hydro}\n")
        stream.write(f"# tau_hyd_fm = {tau_hyd}\n")
        stream.write(f"# tau_min_fm = {tau_min}\n")
        stream.write(f"# eta_over_s = {eta_over_s}\n")
        stream.write(f"# eos_effective_degrees = {QCD_THREE_FLAVOR_EFFECTIVE_DEGREES}\n")
        stream.write(f"# eos_factor = {eos_factor}\n")
        stream.write(f"# attractor = {attractor_label}\n")
        stream.write(
            f"# plan_b_reference = {PLAN_B_REFERENCE}; "
            f"attractor_data_DOI:{PUBLISHED_QCD_ATTRACTOR_DOI}\n"
        )
        stream.write(f"# attractor_C_inf = {QCD_KINETIC_ATTRACTOR_C_INF}\n")
        stream.write(f"# hbarc_GeV_fm = {HBARC_GEV_FM}\n")
        stream.write(f"# viscous_anchor = {int(viscous_anchor)}\n")
        stream.write("# " + "\t".join(columns) + "\n")
        for tau in tau_grid:
            if tau < tau_min or tau >= tau_hyd - 1e-6:
                continue
            for cell in cells:
                temperature_hyd = float(cell["temperature_hyd"])
                temp_eff, omega, e_value = attractor_temperature(
                    tau=tau,
                    tau_hyd=tau_hyd,
                    temperature_hyd=temperature_hyd,
                    eta_over_s=eta_over_s,
                    attractor=attractor,
                    viscous_anchor=viscous_anchor,
                )
                scale = min(max(tau / tau_hyd, 0.0), 1.0)
                vx = float(cell["vx_hyd"]) * scale
                vy = float(cell["vy_hyd"]) * scale
                valid = 1 if temperature_hyd > 0 else 0
                values = [
                    event_id,
                    f"{tau:.8g}",
                    int(cell["ix"]),
                    int(cell["iy"]),
                    f"{float(cell['x_fm']):.8g}",
                    f"{float(cell['y_fm']):.8g}",
                    f"{temp_eff:.8g}",
                    f"{eos_factor * temp_eff ** 4:.8g}",
                    f"{vx:.8g}",
                    f"{vy:.8g}",
                    valid,
                    f"{temperature_hyd:.8g}",
                    f"{float(cell['vx_hyd']):.8g}",
                    f"{float(cell['vy_hyd']):.8g}",
                    f"{omega:.8g}",
                    f"{e_value:.8g}",
                ]
                stream.write("\t".join(str(value) for value in values) + "\n")


def copy_hydro_inputs(hydro_dir: Path, run_dir: Path) -> None:
    for name in ["evolution_all_xyeta.dat", "NcollList.dat"]:
        shutil.copy2(hydro_dir / name, run_dir / name)
    for name in ["README_staged_event.txt", "music_input", "run.log"]:
        if (hydro_dir / name).exists():
            shutil.copy2(hydro_dir / name, run_dir / name)


def run_variant(
    *,
    variant: str,
    run_dir: Path,
    main_bin: Path,
    template: Path,
    hydro_dir: Path | None,
    seed: int,
    events: int,
    aa: bool,
    pthat_min: float,
    pthat_max: float,
    pdf_mode: str,
    lhapdf_set: str,
    timeout_s: int,
    energy_loss_alpha: float,
    broadening_k: float,
    use_prehydro: bool,
    prehydro_tau_grid: list[float],
    prehydro_tau_min: float,
    prehydro_eta_over_s: float,
    prehydro_eos_factor: float,
    prehydro_attractor_label: str,
    prehydro_attractor: Callable[[float], float],
    prehydro_viscous_anchor: bool,
    event_id: int,
    task_id: int,
    hydro_slot: int,
    hydro_ncoll: int,
    hydro_payload_sha256: str,
    centrality: str,
) -> dict[str, str | int | float]:
    run_dir.mkdir(parents=True, exist_ok=True)
    write_pythia_card(template, run_dir / "setup_pythia.cmnd", pthat_min, pthat_max, pdf_mode, lhapdf_set)
    if aa and hydro_dir is not None:
        copy_hydro_inputs(hydro_dir, run_dir)

    prehydro_file = ""
    if use_prehydro:
        prehydro_file = "prehydro_table.tsv"
        generate_planb_prehydro(
            reference_hydro=run_dir / "evolution_all_xyeta.dat",
            output_path=run_dir / prehydro_file,
            event_id=event_id,
            tau_grid=prehydro_tau_grid,
            tau_min=prehydro_tau_min,
            eta_over_s=prehydro_eta_over_s,
            eos_factor=prehydro_eos_factor,
            attractor_label=prehydro_attractor_label,
            attractor=prehydro_attractor,
            viscous_anchor=prehydro_viscous_anchor,
        )

    write_input(
        run_dir / "hybrid_input.dat",
        seed=seed,
        events=events,
        aa=aa,
        energy_loss_alpha=energy_loss_alpha,
        broadening_k=broadening_k,
        use_prehydro=use_prehydro,
        prehydro_file=prehydro_file or "prehydro_table.tsv",
    )

    start = time.time()
    rc = 124
    timed_out = 0
    try:
        rc = run_main(run_dir, main_bin, timeout_s)
    except subprocess.TimeoutExpired:
        timed_out = 1
    elapsed = time.time() - start

    hadron_output = run_dir / "HYBRID_Hadrons.out"
    if rc == 0 and (not hadron_output.exists() or hadron_output.stat().st_size == 0):
        print(f"ERROR: empty or missing {hadron_output}", file=sys.stderr, flush=True)
        rc = 20

    result = {
        "variant": variant,
        "use_prehydro": int(use_prehydro),
        "energy_loss_alpha": energy_loss_alpha,
        "broadening_k": broadening_k,
        "prehydro_file": prehydro_file,
        "returncode": rc,
        "timeout": timed_out,
        "seconds": elapsed,
        "dir": str(run_dir),
    }
    row = (
        "variant\tuse_prehydro\tenergy_loss_alpha\tbroadening_k\tprehydro_file\t"
        "returncode\ttimeout\tseconds\tdir\t"
        "task_id\tseed\tcentrality\thydro_slot\thydro_event_id\thydro_ncoll\t"
        "hydro_payload_sha256\n"
        f"{variant}\t{int(use_prehydro)}\t{energy_loss_alpha}\t"
        f"{broadening_k}\t{prehydro_file}\t{rc}\t{timed_out}\t"
        f"{elapsed:.3f}\t{run_dir}\t{task_id}\t{seed}\t{centrality}\t{hydro_slot}\t"
        f"{event_id}\t{hydro_ncoll}\t{hydro_payload_sha256}\n"
    )
    (run_dir / "summary.tsv").write_text(row)
    return result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", choices=["aa", "pp"], required=True)
    ap.add_argument("--task-id", type=int, required=True)
    ap.add_argument("--seed-offset", type=int, required=True)
    ap.add_argument("--events", type=int, required=True)
    ap.add_argument("--pthat-min", type=float, default=4.0)
    ap.add_argument("--pthat-max", type=float, default=-1.0)
    ap.add_argument("--pdf-mode", choices=["native_npdf", "off", "lhapdf"], default="native_npdf")
    ap.add_argument("--lhapdf-set", default="EPPS21nlo_CT18Anlo_O16/0")
    ap.add_argument("--timeout-s", type=int, default=7200)
    ap.add_argument("--run-name", required=True)
    ap.add_argument("--aa-centrality-index", type=int, default=-1, help="if nonnegative, force every AA task to this staged centrality index")
    ap.add_argument(
        "--aa-task-manifest",
        type=Path,
        help="v2 task-to-hydro assignment manifest (required by v2 AA jobs)",
    )
    prehydro_mode = ap.add_mutually_exclusive_group()
    prehydro_mode.add_argument(
        "--run-prehydro-pair",
        action="store_true",
        help="for AA tasks, run no-prehydro and prehydro variants in the same job",
    )
    prehydro_mode.add_argument(
        "--run-prehydro-only",
        action="store_true",
        help="for AA tasks, run only the Plan-B prehydro variant",
    )
    ap.add_argument(
        "--no-prehydro-alpha",
        type=float,
        default=DEFAULT_NO_PREHYDRO_ALPHA,
        help="HYBRID mode-0 stopping-strength alpha for the no-prehydro baseline",
    )
    ap.add_argument(
        "--prehydro-alpha",
        type=float,
        default=DEFAULT_PREHYDRO_ALPHA,
        help="HYBRID mode-0 stopping-strength alpha for the Plan-B leg",
    )
    ap.add_argument(
        "--additional-prehydro-alpha",
        type=float,
        action="append",
        default=[],
        help=(
            "additional Plan-B alpha to run against the same no-prehydro "
            "baseline; may be repeated"
        ),
    )
    ap.add_argument(
        "--broadening-k",
        type=float,
        default=DEFAULT_BROADENING_K,
        help="HYBRID transverse-broadening kappa shared by both paired legs",
    )
    ap.add_argument("--prehydro-tau-min", type=float, default=0.24)
    ap.add_argument("--prehydro-tau-grid", default=DEFAULT_PREHYDRO_TAU_GRID)
    ap.add_argument(
        "--prehydro-eta-over-s",
        type=float,
        default=PUBLISHED_PREHYDRO_ETA_OVER_S,
    )
    ap.add_argument("--prehydro-eos-factor", type=float, default=QCD_CONFORMAL_EOS_FACTOR)
    ap.add_argument(
        "--prehydro-attractor-table",
        type=Path,
        default=None,
        help="published two-column omega,E(omega) table; required for either prehydro mode",
    )
    ap.add_argument("--prehydro-viscous-anchor", dest="prehydro_viscous_anchor", action="store_true", default=True)
    ap.add_argument("--no-prehydro-viscous-anchor", dest="prehydro_viscous_anchor", action="store_false")
    args = ap.parse_args()

    for name, value in (
        ("--no-prehydro-alpha", args.no_prehydro_alpha),
        ("--prehydro-alpha", args.prehydro_alpha),
        ("--broadening-k", args.broadening_k),
        *(
            ("--additional-prehydro-alpha", value)
            for value in args.additional_prehydro_alpha
        ),
    ):
        if not math.isfinite(value) or value < 0.0:
            ap.error(f"{name} must be finite and nonnegative, got {value}")

    work = Path.cwd()
    main_bin = Path(os.environ.get("MMLI_BIN", work / "bin" / "main"))
    template = Path(os.environ.get("MMLI_PYTHIA_TEMPLATE", work / "runtime" / "setup_pythia.cmnd"))
    hydro_root = Path(os.environ.get("OO_HYDRO_ROOT", work / "runtime" / "staged_hydro"))
    prehydro_tau_grid = parse_tau_grid(args.prehydro_tau_grid)
    if not math.isclose(args.prehydro_eos_factor, QCD_CONFORMAL_EOS_FACTOR, rel_tol=0.0, abs_tol=1e-12):
        ap.error(
            "--prehydro-eos-factor must use the three-flavor conformal value "
            f"{QCD_CONFORMAL_EOS_FACTOR:.15g}"
        )
    prehydro_enabled = args.run_prehydro_pair or args.run_prehydro_only
    if args.kind != "aa" and prehydro_enabled:
        ap.error("prehydro modes are only valid for AA tasks")
    if args.kind == "aa" and prehydro_enabled:
        if not math.isclose(
            args.prehydro_eta_over_s,
            PUBLISHED_PREHYDRO_ETA_OVER_S,
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            ap.error(
                "--prehydro-eta-over-s must be 0.12 for the published "
                "arXiv:2509.19430v2 Plan B prescription"
            )
        if not args.prehydro_viscous_anchor:
            ap.error(
                "--no-prehydro-viscous-anchor is incompatible with the "
                "published arXiv:2509.19430v2 Plan B prescription"
            )
        if args.prehydro_attractor_table is None:
            ap.error("--prehydro-attractor-table is required for prehydro production")
        prehydro_attractor_label, prehydro_attractor = load_attractor_table(
            args.prehydro_attractor_table
        )
    else:
        prehydro_attractor_label = "disabled"

        def prehydro_attractor(_: float) -> float:
            raise RuntimeError("prehydro attractor called for a non-prehydro run")

    seed = args.seed_offset + args.task_id
    if args.kind == "aa":
        if args.aa_task_manifest is not None:
            assignment = load_aa_task_assignment(args.aa_task_manifest, args.task_id)
            if int(assignment["hard_seed"]) != seed:
                raise ValueError(
                    f"task {args.task_id}: manifest hard seed {assignment['hard_seed']} "
                    f"does not match seed offset result {seed}"
                )
            hydro_index = int(assignment["hydro_slot"])
            cent = "C0-5"
            hydro_dir = hydro_root / str(assignment["hydro_dir"])
            event_id = int(assignment["hydro_event_id"])
            hydro_ncoll = int(assignment["hydro_ncoll"])
            hydro_payload_sha256 = str(assignment["hydro_payload_sha256"])
            staged_metadata = read_staged_hydro_metadata(hydro_dir)
            expected_metadata = {
                "hydro_slot": hydro_index,
                "event_id": event_id,
                "ncoll_positions": hydro_ncoll,
            }
            if staged_metadata != expected_metadata:
                raise ValueError(
                    f"task {args.task_id}: staged hydro metadata {staged_metadata} "
                    f"does not match manifest {expected_metadata}"
                )
            run_hydro_label = f"hydro_{hydro_index:03d}_{cent}_event_{event_id:05d}"
        elif args.aa_centrality_index >= 0:
            if args.aa_centrality_index >= len(CENTRALITIES):
                raise ValueError(f"--aa-centrality-index must be 0..{len(CENTRALITIES) - 1}")
            hydro_index = args.aa_centrality_index
            cent = CENTRALITIES[hydro_index]
            hydro_dir = hydro_root / f"{cent}_idx0"
            event_id = event_id_from_readme(hydro_dir, args.task_id)
            hydro_ncoll = -1
            hydro_payload_sha256 = ""
            run_hydro_label = f"hydro_{hydro_index:02d}_{cent}"
        else:
            hydro_index = args.task_id % len(CENTRALITIES)
            cent = CENTRALITIES[hydro_index]
            hydro_dir = hydro_root / f"{cent}_idx0"
            event_id = event_id_from_readme(hydro_dir, args.task_id)
            hydro_ncoll = -1
            hydro_payload_sha256 = ""
            run_hydro_label = f"hydro_{hydro_index:02d}_{cent}"
        if not (hydro_dir / "evolution_all_xyeta.dat").exists():
            raise FileNotFoundError(hydro_dir / "evolution_all_xyeta.dat")
        base_run_dir = work / args.run_name / "aa" / run_hydro_label / f"task_{args.task_id:05d}"
        variants = build_aa_variants(
            base_run_dir=base_run_dir,
            task_id=args.task_id,
            run_prehydro_pair=args.run_prehydro_pair,
            no_prehydro_alpha=args.no_prehydro_alpha,
            prehydro_alpha=args.prehydro_alpha,
            additional_prehydro_alphas=args.additional_prehydro_alpha,
            run_prehydro_only=args.run_prehydro_only,
        )
    else:
        hydro_index = -1
        cent = "pp_reference"
        hydro_dir = None
        event_id = args.task_id
        hydro_ncoll = -1
        hydro_payload_sha256 = ""
        variants = [
            (
                "pp_reference",
                work / args.run_name / "pp" / f"task_{args.task_id:05d}",
                False,
                args.no_prehydro_alpha,
            )
        ]

    results = []
    for variant, run_dir, use_prehydro, energy_loss_alpha in variants:
        results.append(
            run_variant(
                variant=variant,
                run_dir=run_dir,
                main_bin=main_bin,
                template=template,
                hydro_dir=hydro_dir,
                seed=seed,
                events=args.events,
                aa=args.kind == "aa",
                pthat_min=args.pthat_min,
                pthat_max=args.pthat_max,
                pdf_mode=args.pdf_mode,
                lhapdf_set=args.lhapdf_set,
                timeout_s=args.timeout_s,
                energy_loss_alpha=energy_loss_alpha,
                broadening_k=args.broadening_k,
                use_prehydro=use_prehydro,
                prehydro_tau_grid=prehydro_tau_grid,
                prehydro_tau_min=args.prehydro_tau_min,
                prehydro_eta_over_s=args.prehydro_eta_over_s,
                prehydro_eos_factor=args.prehydro_eos_factor,
                prehydro_attractor_label=prehydro_attractor_label,
                prehydro_attractor=prehydro_attractor,
                prehydro_viscous_anchor=args.prehydro_viscous_anchor,
                event_id=event_id,
                task_id=args.task_id,
                hydro_slot=hydro_index,
                hydro_ncoll=hydro_ncoll,
                hydro_payload_sha256=hydro_payload_sha256,
                centrality=cent,
            )
        )

    header = (
        "kind\ttask_id\tseed\tevents\tcentrality\thydro_index\thydro_event_id\t"
        "hydro_ncoll\thydro_payload_sha256\tvariant\tuse_prehydro\t"
        "energy_loss_alpha\tbroadening_k\tprehydro_file\treturncode\ttimeout\t"
        "seconds\tdir\n"
    )
    rows = []
    for result in results:
        rows.append(
            f"{args.kind}\t{args.task_id}\t{seed}\t{args.events}\t{cent}\t{hydro_index}\t"
            f"{event_id}\t{hydro_ncoll}\t{hydro_payload_sha256}\t"
            f"{result['variant']}\t{result['use_prehydro']}\t"
            f"{result['energy_loss_alpha']}\t{result['broadening_k']}\t"
            f"{result['prehydro_file']}\t"
            f"{result['returncode']}\t{result['timeout']}\t{float(result['seconds']):.3f}\t{result['dir']}"
        )
    text = header + "\n".join(rows) + "\n"
    summary_dir = Path(str(results[0]["dir"])).parent
    summary_kind = "prehydro_only" if args.run_prehydro_only else "pair"
    (summary_dir / f"task_{args.task_id:05d}_{summary_kind}_summary.tsv").write_text(text)
    print(text, end="")
    return 0 if all(int(r["returncode"]) == 0 and int(r["timeout"]) == 0 for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
