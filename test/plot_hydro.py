#!/usr/bin/env python3
"""
Plot temperature profile T(x, y) at several tau slices from
hydroinfoPlaintxtHuichaoFormat.dat
Columns: x  y  tau  e  T  vx  vy
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

HYDRO_FILE = "./hydroinfoPlaintxtHuichaoFormat.dat"

# Temperature conversion: raw file is in fm^-1 (natural units)
# T[GeV] = T[fm^-1] * hbar*c,  hbar*c = 0.197327 GeV·fm
HBARC = 0.197327

# ── 1. Stream through file, collect unique tau values and build a dict ─────────
print("Reading hydro file …", flush=True)

# First pass: collect unique taus
unique_taus = set()
with open(HYDRO_FILE) as f:
    for line in f:
        vals = line.split()
        if len(vals) >= 5:
            unique_taus.add(float(vals[2]))

sorted_taus = sorted(unique_taus)
print(f"  Found {len(sorted_taus)} unique tau values: {sorted_taus[0]:.2f} … {sorted_taus[-1]:.2f} fm")

# Choose ~6 representative tau snapshots (evenly spaced in the list)
n_panels = 6
idx_choices = np.linspace(0, len(sorted_taus) - 1, n_panels, dtype=int)
chosen_taus = [sorted_taus[i] for i in idx_choices]
chosen_set = set(chosen_taus)
print(f"  Plotting taus: {[f'{t:.2f}' for t in chosen_taus]}")

# Second pass: collect (x, y, T) for chosen taus only
data = {tau: [] for tau in chosen_taus}
with open(HYDRO_FILE) as f:
    for line in f:
        vals = line.split()
        if len(vals) < 5:
            continue
        x, y, tau = float(vals[0]), float(vals[1]), float(vals[2])
        if tau in chosen_set:
            T = float(vals[4]) * HBARC
            data[tau].append((x, y, T))

print("  Data loaded.", flush=True)

# ── 2. Build 2-D grids and plot ────────────────────────────────────────────────
fig = plt.figure(figsize=(18, 10))
fig.suptitle("Temperature profile T(x, y) [GeV] at different proper times τ",
             fontsize=14, fontweight="bold")

gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.4, wspace=0.35)

for panel_idx, tau in enumerate(chosen_taus):
    arr = np.array(data[tau])
    xs = np.unique(arr[:, 0])
    ys = np.unique(arr[:, 1])

    # Build 2-D temperature grid (rows = y, cols = x)
    T_grid = np.zeros((len(ys), len(xs)))
    xi_map = {v: i for i, v in enumerate(xs)}
    yi_map = {v: i for i, v in enumerate(ys)}
    for x, y, T in arr:
        T_grid[yi_map[y], xi_map[x]] = T

    ax = fig.add_subplot(gs[panel_idx // 3, panel_idx % 3])
    im = ax.pcolormesh(xs, ys, T_grid, cmap="inferno", shading="auto")
    plt.colorbar(im, ax=ax, label="T [GeV]")
    ax.set_title(f"τ = {tau:.2f} fm", fontsize=11)
    ax.set_xlabel("x [fm]")
    ax.set_ylabel("y [fm]")
    ax.set_aspect("equal")

out_file = "hydro_temperature_profile.pdf"
fig.savefig(out_file, bbox_inches="tight")
print(f"Saved: {out_file}")

out_png = "hydro_temperature_profile.png"
fig.savefig(out_png, dpi=150, bbox_inches="tight")
print(f"Saved: {out_png}")
plt.show()
