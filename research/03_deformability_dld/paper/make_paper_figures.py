"""
Generate the extra figures requested for the M&N paper:

  Figure 1 — snapshots: initial state (top) and final state (bottom)
             of the peak-signal cell, with pillar circles overlaid.
  Figure 4 — trajectory overlay: soft (purple) vs stiff (cyan) trails
             across the full pillar array, peak-signal cell.

Reads sweep_out/G0.030_S0.400/history.npz and writes to
paper/figures/.

Usage:
    python make_paper_figures.py
"""
from __future__ import annotations

import pathlib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


# ------------------------------------------------------------------
# Geometry — must match sweep.py exactly
# ------------------------------------------------------------------
NX, NY = 400, 80
N_PILLAR_COLS  = 8
N_PILLAR_ROWS  = 4
PILLAR_R       = 4.0
PILLAR_DX      = 35.0
PILLAR_DY      = 14.0
ROW_SHIFT      = 3.5
ARRAY_X_START  = 50.0
ARRAY_Y_CENTRE = NY / 2.0
R_PART         = 3.0

# Visual style
SOFT_COLOR  = "#8e44ad"    # purple
STIFF_COLOR = "#16a085"    # teal-cyan
PILLAR_COLOR = "#555555"

PEAK_CELL = pathlib.Path(__file__).resolve().parents[1] / "sweep_out" / "G0.030_S0.400"


# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------
def pillar_positions():
    """Return list of (cx, cy) for all pillars in the array (matches
    the geometry built by run.py / sweep.py)."""
    centres = []
    for col in range(N_PILLAR_COLS):
        for row in range(N_PILLAR_ROWS):
            cx = ARRAY_X_START + col * PILLAR_DX
            cy = (ARRAY_Y_CENTRE
                  + (row - (N_PILLAR_ROWS - 1) / 2.0) * PILLAR_DY
                  + col * ROW_SHIFT)
            if cy < 6 or cy > NY - 6:
                continue
            centres.append((cx, cy))
    return centres


def draw_pillars(ax):
    for cx, cy in pillar_positions():
        ax.add_patch(Circle((cx, cy), PILLAR_R,
                            facecolor=PILLAR_COLOR, edgecolor="black",
                            linewidth=0.4, alpha=0.85, zorder=2))


def draw_walls(ax):
    """Highlight the y=0 and y=NY walls."""
    ax.axhline(0.5, color="black", linewidth=1.5, zorder=1)
    ax.axhline(NY - 0.5, color="black", linewidth=1.5, zorder=1)


def style_axes(ax, title=""):
    ax.set_xlim(0, NX)
    ax.set_ylim(0, NY)
    ax.set_aspect("equal")
    ax.set_xlabel("x  (lattice units)")
    ax.set_ylabel("y  (lattice units)")
    ax.set_title(title)
    ax.grid(False)


# ------------------------------------------------------------------
# Figure 1 — initial vs final state
# ------------------------------------------------------------------
def figure_initial_final(history_path, out_path):
    h = np.load(history_path, allow_pickle=True)
    traj  = h["traj_xy"]            # (n_frames, n_particles, 2)
    types = h["types"]              # (n_particles,)
    soft_mask  = (types == 0)
    stiff_mask = (types == 1)

    pos0 = traj[0]
    posN = traj[-1]
    # Periodic-x wrap correction for visualisation: keep particles in [0, NX)
    posN_vis = posN.copy()
    posN_vis[:, 0] = posN_vis[:, 0] % NX

    fig, axes = plt.subplots(2, 1, figsize=(11, 6.2))

    # --- Top: t = 0 -------------------------------------------------
    ax = axes[0]
    draw_pillars(ax)
    draw_walls(ax)
    for x, y in pos0[soft_mask]:
        ax.add_patch(Circle((x, y), R_PART, facecolor=SOFT_COLOR,
                             edgecolor="black", linewidth=0.4,
                             alpha=0.85, zorder=3))
    for x, y in pos0[stiff_mask]:
        ax.add_patch(Circle((x, y), R_PART, facecolor=STIFF_COLOR,
                             edgecolor="black", linewidth=0.4,
                             alpha=0.85, zorder=3))
    style_axes(ax, title="(a) Initial state, $t = 0$")

    # Legend (manual)
    ax.plot([], [], "o", color=SOFT_COLOR,  markersize=10,
            markeredgecolor="black", markeredgewidth=0.5,
            label=fr"soft ($G_s = 0.030$, $n=$ {soft_mask.sum()})")
    ax.plot([], [], "o", color=STIFF_COLOR, markersize=10,
            markeredgecolor="black", markeredgewidth=0.5,
            label=fr"stiff ($G_s = 0.400$, $n=$ {stiff_mask.sum()})")
    ax.plot([], [], "o", color=PILLAR_COLOR, markersize=10,
            label="pillar")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.92)

    # --- Bottom: t = T (final) --------------------------------------
    ax = axes[1]
    draw_pillars(ax)
    draw_walls(ax)
    for x, y in posN_vis[soft_mask]:
        ax.add_patch(Circle((x, y), R_PART, facecolor=SOFT_COLOR,
                             edgecolor="black", linewidth=0.4,
                             alpha=0.85, zorder=3))
    for x, y in posN_vis[stiff_mask]:
        ax.add_patch(Circle((x, y), R_PART, facecolor=STIFF_COLOR,
                             edgecolor="black", linewidth=0.4,
                             alpha=0.85, zorder=3))
    style_axes(ax, title=f"(b) Final state, $t = 20\,000$  (after {traj.shape[0] - 1} snapshots)")

    fig.suptitle(
        r"One representative cell:  $G_s^{\rm soft} = 0.030$, "
        r"$G_s^{\rm stiff} = 0.400$, $\mathrm{Ca}_{\rm ratio} = 13.3$",
        fontsize=10, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


# ------------------------------------------------------------------
# Figure 4 — trajectory overlay (soft vs stiff trails)
# ------------------------------------------------------------------
def figure_trajectories(history_path, out_path):
    h = np.load(history_path, allow_pickle=True)
    traj  = h["traj_xy"]
    types = h["types"]
    soft_mask  = (types == 0)
    stiff_mask = (types == 1)

    # Periodic-x un-wrap for clean trails: convert wrapped jumps into
    # continuous cumulative-x trajectories
    traj_unwrap = traj.copy()
    for p in range(traj.shape[1]):
        x = traj_unwrap[:, p, 0]
        dx = np.diff(x)
        jumps = np.where(dx < -NX / 2)[0]   # wrap-around
        for j in jumps:
            traj_unwrap[j + 1:, p, 0] += NX

    fig, ax = plt.subplots(figsize=(11, 4.4))
    draw_pillars(ax)
    draw_walls(ax)

    # Plot trails — semi-transparent thin lines from initial to final
    for p in np.where(soft_mask)[0]:
        ax.plot(traj_unwrap[:, p, 0] % NX, traj_unwrap[:, p, 1],
                color=SOFT_COLOR, linewidth=0.5, alpha=0.5, zorder=4)
    for p in np.where(stiff_mask)[0]:
        ax.plot(traj_unwrap[:, p, 0] % NX, traj_unwrap[:, p, 1],
                color=STIFF_COLOR, linewidth=0.5, alpha=0.5, zorder=4)

    # Initial markers
    ax.scatter(traj[0, soft_mask, 0], traj[0, soft_mask, 1],
                marker="o", s=12, color=SOFT_COLOR,
                edgecolor="black", linewidths=0.3, zorder=5,
                label=fr"soft initial ($n=$ {soft_mask.sum()})")
    ax.scatter(traj[0, stiff_mask, 0], traj[0, stiff_mask, 1],
                marker="s", s=12, color=STIFF_COLOR,
                edgecolor="black", linewidths=0.3, zorder=5,
                label=fr"stiff initial ($n=$ {stiff_mask.sum()})")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.92)

    style_axes(ax, title=r"Trajectory overlay: peak-signal cell"
                          r" ($G_s^{\rm soft} = 0.030$, $G_s^{\rm stiff} = 0.400$)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def main():
    history = PEAK_CELL / "history.npz"
    if not history.exists():
        raise SystemExit(f"history.npz not found at {history}")
    out_dir = pathlib.Path(__file__).resolve().parent / "figures"
    out_dir.mkdir(exist_ok=True)
    figure_initial_final(history, out_dir / "fig_snapshots.png")
    figure_trajectories(history, out_dir / "fig_trajectories.png")


if __name__ == "__main__":
    main()
