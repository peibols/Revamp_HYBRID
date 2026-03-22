#!/usr/bin/env python3
"""
Plot velocity vector field (vx, vy) + temperature contour at several tau
snapshots from hydroinfoPlaintxtHuichaoFormat.dat
Columns: x  y  tau  e  T  vx  vy

Arrow tips show flow direction; arrow length ∝ |v|/v_max, scaled so the
longest arrows span one quiver cell (STEP × dx).
Temperature contour at T = Tc = 0.170 GeV marks the phase boundary.
"""

import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

HYDRO_FILE = "./hydroinfoPlaintxtHuichaoFormat.dat"
HBARC = 0.197327  # GeV·fm: converts fm^-1 -> GeV
TC = 0.145        # GeV – deconfinement temperature contour
STEP = 5          # quiver downsampling: 1 arrow per STEP grid points
DX   = 0.3        # fm – grid spacing
ZOOM = 14.0       # fm – half-width of plotted region

# ── 1. Collect unique tau values ───────────────────────────────────────────────
print("Scanning tau values …", flush=True)
unique_taus = set()
with open(HYDRO_FILE) as f:
    for line in f:
        cols = line.split()
        if len(cols) >= 7:
            unique_taus.add(float(cols[2]))

sorted_taus = sorted(unique_taus)
print(f"  {len(sorted_taus)} tau slices: {sorted_taus[0]:.2f} – {sorted_taus[-1]:.2f} fm")

# Choose 6 representative snapshots (skip tau=0.6 since v=0 everywhere)
n_panels = 6
chosen_taus = [sorted_taus[i] for i in np.linspace(1, len(sorted_taus)-1, n_panels, dtype=int)]
chosen_set  = set(chosen_taus)
print(f"  Plotting taus: {[f'{t:.2f}' for t in chosen_taus]}")

# ── 2. Load data ───────────────────────────────────────────────────────────────
data = {tau: [] for tau in chosen_taus}
with open(HYDRO_FILE) as f:
    for line in f:
        cols = line.split()
        if len(cols) < 7:
            continue
        tau = float(cols[2])
        if tau in chosen_set:
            data[tau].append([float(cols[0]), float(cols[1]),
                              float(cols[4]) * HBARC,       # T [GeV]
                              float(cols[5]), float(cols[6])])  # vx, vy
print("  Data loaded.", flush=True)

# ── 3. Build grids and plot ────────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 11))
fig.suptitle(
    "Transverse flow velocity |v| (colour) + direction (arrows)\n"
    f"Dashed contour: T = {TC:.3f} GeV (phase boundary)",
    fontsize=13, fontweight="bold")

gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.50, wspace=0.30)

for panel_idx, tau in enumerate(chosen_taus):
    arr = np.array(data[tau])   # columns: x, y, T, vx, vy

    xs = np.unique(arr[:,0])
    ys = np.unique(arr[:,1])
    nx, ny = len(xs), len(ys)

    xi_map = {v: i for i, v in enumerate(xs)}
    yi_map = {v: i for i, v in enumerate(ys)}

    T_grid  = np.zeros((ny, nx))
    Vx_grid = np.zeros((ny, nx))
    Vy_grid = np.zeros((ny, nx))

    for row in arr:
        x, y, T, vx, vy = row
        i, j = yi_map[y], xi_map[x]
        T_grid[i,j]  = T
        Vx_grid[i,j] = vx
        Vy_grid[i,j] = vy

    Vmag = np.sqrt(Vx_grid**2 + Vy_grid**2)
    vmax = Vmag.max()

    ax = fig.add_subplot(gs[panel_idx // 3, panel_idx % 3])

    # ── background: |v| ───────────────────────────────────────────────────────
    im = ax.pcolormesh(xs, ys, Vmag, cmap="plasma",
                       vmin=0, vmax=min(vmax, 0.99), shading="auto")
    cb = plt.colorbar(im, ax=ax, label="|v| / c", pad=0.02)
    cb.ax.tick_params(labelsize=8)

    # ── temperature contour at Tc ─────────────────────────────────────────────
    ax.contour(xs, ys, T_grid, levels=[TC],
               colors="cyan", linewidths=1.2, linestyles="--")

    # ── velocity arrows ────────────────────────────────────────────────────────
    X2d, Y2d = np.meshgrid(xs, ys)
    # Normalise to unit vectors so arrow lengths are all equal (pure direction)
    # but zero-velocity cells get no arrow via a mask
    Vmag_ds   = Vmag[::STEP, ::STEP]
    mask      = Vmag_ds > 0.01          # suppress arrows where fluid is at rest
    Vx_ds = np.where(mask, Vx_grid[::STEP, ::STEP] / (Vmag_ds + 1e-12), 0.0)
    Vy_ds = np.where(mask, Vy_grid[::STEP, ::STEP] / (Vmag_ds + 1e-12), 0.0)

    # scale: 1 normalised unit → STEP*DX fm  (one arrow spacing per max arrow)
    arrow_scale = 1.0 / (STEP * DX)    # data units per normalised unit
    ax.quiver(X2d[::STEP, ::STEP], Y2d[::STEP, ::STEP],
              Vx_ds, Vy_ds,
              color="white", alpha=0.80,
              scale=arrow_scale, scale_units="xy",
              width=0.005, headwidth=3.5, headlength=4)

    ax.set_xlim(-ZOOM, ZOOM)
    ax.set_ylim(-ZOOM, ZOOM)
    ax.set_title(f"τ = {tau:.2f} fm   (v$_{{\\rm max}}$={vmax:.2f}c)", fontsize=10)
    ax.set_xlabel("x [fm]", fontsize=9)
    ax.set_ylabel("y [fm]", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_aspect("equal")

for ext in ("pdf", "png"):
    out = f"hydro_velocity_profile.{ext}"
    fig.savefig(out, bbox_inches="tight", dpi=(150 if ext == "png" else None))
    print(f"Saved: {out}")

plt.close("all")
