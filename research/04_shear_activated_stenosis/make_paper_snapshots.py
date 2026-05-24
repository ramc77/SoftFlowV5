"""
Generate publication-grade snapshot figures from a single cell's output.
Default: sev50_th2 (sweet-spot cell).

Reads:
  - sweep_out/<cell>/history.npz                  (drug-delivery time series)
  - sweep_out/<cell>/vtk_stenosis/particle_data.csv (per-step carrier positions)
  - sweep_out/<cell>/vtk_stenosis/global_stats.csv  (fluid + drug aggregate stats)

Writes:
  - fig_schematic.png       — annotated geometry diagram
  - fig_snapshots.png       — 3-panel spatial state at t = 0, 100k, 200k
  - fig_trajectory.png      — time series of deposition / release / off-target
"""
import sys
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

# Configuration — match run.py exactly
NX, NY            = 400, 80
STENOSIS_R        = 50.0
STENOSIS_X        = 200
TARGET_X_LO       = 140
TARGET_X_HI       = 280
OFFTARGET_X_UP    = (10, 110)
OFFTARGET_X_DN    = (310, 390)
OFFTARGET_BAND    = 3
R_CARRIER         = 2.0

# Defaults for the sweet-spot cell
DEFAULT_CELL      = "sev50_th2"
DEFAULT_SEVERITY  = 50
DEFAULT_GAMMA_TH  = 1e-3

# Colours
COL_CARRIER       = "#c0392b"
COL_OBSTACLE      = "#7f8c8d"
COL_TARGET        = "#16a085"
COL_OFFTARGET     = "#3498db"


def throat_extent(severity_pct):
    if severity_pct == 0:
        return (0, NY)
    p = int((severity_pct / 100.0) * (NY / 2.0))
    return (p, NY - p)


def obstacle_centres(severity_pct):
    if severity_pct == 0:
        return []
    p = (severity_pct / 100.0) * (NY / 2.0)
    return [(STENOSIS_X, -STENOSIS_R + p), (STENOSIS_X, NY + STENOSIS_R - p)]


def draw_geometry(ax, severity, alpha_target=0.18, alpha_off=0.13, draw_zones=True):
    """Render the channel, stenosis obstacles, target + off-target zones."""
    # Channel walls (top + bottom)
    ax.axhline(0, color="black", linewidth=1.5)
    ax.axhline(NY, color="black", linewidth=1.5)
    # Stenosis obstacles
    for (cx, cy) in obstacle_centres(severity):
        ax.add_patch(Circle((cx, cy), STENOSIS_R, color=COL_OBSTACLE,
                             zorder=2))
    # Target zone (throat-region absorber)
    if draw_zones and severity > 0:
        j_lo, j_hi = throat_extent(severity)
        ax.add_patch(Rectangle(
            (TARGET_X_LO, j_lo), TARGET_X_HI - TARGET_X_LO, j_hi - j_lo,
            facecolor=COL_TARGET, alpha=alpha_target, edgecolor=COL_TARGET,
            linewidth=1.2, linestyle="--", zorder=1,
            label="target absorber"))
    # Off-target wall zones
    if draw_zones:
        for (x_lo, x_hi) in (OFFTARGET_X_UP, OFFTARGET_X_DN):
            for (y_lo, y_hi) in ((0, OFFTARGET_BAND),
                                  (NY - OFFTARGET_BAND, NY)):
                ax.add_patch(Rectangle(
                    (x_lo, y_lo), x_hi - x_lo, y_hi - y_lo,
                    facecolor=COL_OFFTARGET, alpha=alpha_off,
                    edgecolor=COL_OFFTARGET, linewidth=1.0,
                    zorder=1))
    ax.set_xlim(0, NX)
    ax.set_ylim(-2, NY + 2)
    ax.set_aspect("equal")
    ax.set_xlabel("x (LU)")
    ax.set_ylabel("y (LU)")


def fig_schematic(out_path, severity=DEFAULT_SEVERITY, gamma_th=DEFAULT_GAMMA_TH):
    """Annotated geometry diagram."""
    fig, ax = plt.subplots(figsize=(11, 3.5))
    draw_geometry(ax, severity)
    # Carrier seed region (upstream)
    ax.add_patch(Rectangle((10, 12), 70, NY - 24,
                            facecolor=COL_CARRIER, alpha=0.12,
                            edgecolor=COL_CARRIER, linewidth=1.0,
                            linestyle=":", label="carrier seed zone"))
    # Annotations
    ax.annotate("flow", xy=(370, NY/2), xytext=(320, NY/2),
                arrowprops=dict(arrowstyle="->", lw=2, color="black"),
                fontsize=11, ha="center")
    ax.annotate("stenosis throat", xy=(STENOSIS_X, NY/2),
                xytext=(STENOSIS_X, NY + 12),
                arrowprops=dict(arrowstyle="->", lw=1.0),
                fontsize=10, ha="center")
    ax.annotate("target zone\n(diseased wall, $k_{\\rm target}=0.5$)",
                xy=((TARGET_X_LO + TARGET_X_HI) / 2, NY + 12),
                fontsize=8, color=COL_TARGET, ha="center")
    ax.annotate("off-target\n($k_{\\rm off}=0.05$)",
                xy=(60, -7), fontsize=8, color=COL_OFFTARGET, ha="center")
    ax.set_title(f"System schematic — $\\sigma$ = {severity}%, "
                  f"$\\dot\\gamma_{{\\rm th}}$ = {gamma_th:.0e} LU")
    ax.set_ylim(-15, NY + 25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def fig_snapshots(cell_dir, out_path, severity=DEFAULT_SEVERITY):
    """3-panel spatial state at t = 0, ~100k, 200k."""
    csv = cell_dir / "vtk_stenosis" / "particle_data.csv"
    df = pd.read_csv(csv)
    # Pick three timesteps roughly equally spaced
    all_t = sorted(df["timestep"].unique())
    chosen = [all_t[0], all_t[len(all_t) // 2], all_t[-1]]
    labels = [f"$t$ = {int(t)}" for t in chosen]

    fig, axes = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True)
    for ax, t, label in zip(axes, chosen, labels):
        draw_geometry(ax, severity, alpha_target=0.10, alpha_off=0.08)
        sub = df[df["timestep"] == t]
        for _, row in sub.iterrows():
            ax.add_patch(Circle((row["cx"], row["cy"]), R_CARRIER,
                                 facecolor=COL_CARRIER, edgecolor="black",
                                 linewidth=0.3, alpha=0.85, zorder=3))
        ax.text(0.02, 0.92, label, transform=ax.transAxes,
                fontsize=11, fontweight="bold",
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none"))
    axes[0].set_title(f"Spatial trajectory — sweet-spot cell "
                       f"($\\sigma$ = {severity}%, $\\dot\\gamma_{{\\rm th}}$ = 1e-3 LU)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def fig_trajectory(cell_dir, out_path):
    """Time series of deposition / release / off-target / drug in fluid."""
    h = np.load(cell_dir / "history.npz", allow_pickle=True)
    hist = h["history"]
    t   = hist["time"]
    rel = hist["total_released"]
    tgt = hist["target_absorbed"]
    off = hist["off_target_absorbed"]
    in_fluid = rel - tgt - off  # mass conservation (approximately)
    M_total = 18.0  # 18 carriers × 1.0 each

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    # Top: absolute amounts
    ax1.plot(t, rel, color="black", linewidth=2, label="total released")
    ax1.plot(t, tgt, color=COL_TARGET, linewidth=2,
              label="absorbed at target (stenosis)")
    ax1.plot(t, off, color=COL_OFFTARGET, linewidth=2,
              label="absorbed off-target (healthy wall)")
    ax1.plot(t, in_fluid, color="#9b59b6", linewidth=1.5,
              linestyle="--", label="remaining in fluid")
    ax1.set_ylabel("drug mass (LU units)")
    ax1.set_title(f"Single-cell deposition trajectory ($\\sigma$ = {h['severity']}%, "
                   f"$\\dot\\gamma_{{\\rm th}}$ = {h['gamma_th']:.0e} LU, "
                   f"200 000 steps ≈ 10 cycles)")
    ax1.legend(loc="center right", fontsize=9, frameon=True)
    ax1.grid(alpha=0.3)

    # Bottom: η_deposit(t) and η_off(t)
    with np.errstate(divide="ignore", invalid="ignore"):
        eta_dep = np.where(rel > 1e-12, tgt / rel, 0.0)
        eta_off = np.where(rel > 1e-12, off / rel, 0.0)
    ax2.plot(t, eta_dep, color=COL_TARGET, linewidth=2,
              label=r"$\eta_{\rm deposit}(t)$")
    ax2.plot(t, eta_off, color=COL_OFFTARGET, linewidth=2,
              label=r"$\eta_{\rm off}(t)$")
    ax2.axhline(0.5, color="gray", linewidth=0.7, linestyle=":")
    ax2.set_ylim(0, 1)
    ax2.set_xlabel("time (LU)")
    ax2.set_ylabel("efficiency (fraction of released)")
    ax2.legend(loc="center right", fontsize=9, frameon=True)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main():
    base = pathlib.Path(__file__).resolve().parent
    cell_dir = base / "sweep_out" / DEFAULT_CELL
    if not cell_dir.exists():
        sys.exit(f"Cell directory {cell_dir} not found")

    fig_schematic(base / "fig_schematic.png")
    fig_snapshots(cell_dir, base / "fig_snapshots.png")
    fig_trajectory(cell_dir, base / "fig_trajectory.png")


if __name__ == "__main__":
    main()
