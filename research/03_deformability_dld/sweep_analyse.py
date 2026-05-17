"""Plot the Project-3 stiffness sweep results.

Reads ``sweep_results.npz`` (produced by ``sweep.py``) and emits
publication-figure candidates:

  fig_sweep_dtheta.png       Δθ heatmap on the (G_s_soft, G_s_stiff)
                              grid — the headline.
  fig_sweep_lane_order.png   Lane-order parameter Φ_lane on the
                              same grid.
  fig_sweep_dtheta_vs_Ca.png Δθ vs Ca-ratio (collapse onto a single
                              dimensionless axis — proof the effect
                              is deformability-controlled).
  fig_sweep_summary.txt      Plain-text table of every cell, easy
                              to drop into a paper's supplement.

If matplotlib isn't installed, only the text summary is written.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np


HERE = pathlib.Path(__file__).resolve().parent


def _grid_from_pairs(soft_vals: np.ndarray, stiff_vals: np.ndarray,
                      d_theta: np.ndarray,
                      sweep_axis: np.ndarray) -> np.ndarray:
    """Re-fold the linear list of cells back into a (N, N) matrix
    indexed by the sweep-axis values. Missing cells stay NaN."""
    n = sweep_axis.size
    grid = np.full((n, n), np.nan, dtype=np.float64)
    for s, t, v in zip(soft_vals, stiff_vals, d_theta):
        i = int(np.argmin(np.abs(sweep_axis - s)))
        j = int(np.argmin(np.abs(sweep_axis - t)))
        grid[i, j] = v
    return grid


def diagonal_baseline(results: dict) -> dict[float, float]:
    """Return G_s_soft → Δθ(soft, soft) baseline for diagonal subtraction.

    The diagonal cell (G_s_soft = G_s_stiff) is a no-contrast control:
    both species are physically identical, so any non-zero Δθ on the
    diagonal is a systematic placement bias (e.g. seeded in different
    inlet sub-regions). We subtract this row-wise baseline from every
    cell to recover the *deformability-only* signal.
    """
    base: dict[float, float] = {}
    for k in range(int(results["n_cells"])):
        s = float(results["G_s_soft"][k])
        t = float(results["G_s_stiff"][k])
        if abs(s - t) < 1e-9:
            base[s] = float(results["d_theta"][k])
    return base


def write_text_summary(out_path: pathlib.Path, results: dict) -> None:
    """Plain-text table — copy-paste into the paper's supplement."""
    base = diagonal_baseline(results)
    lines = []
    lines.append("Project 3 — DLD deformability sweep summary")
    lines.append("=" * 90)
    lines.append(
        f"{'G_s_soft':>10s} {'G_s_stiff':>10s} {'Ca_ratio':>9s} "
        f"{'theta_soft':>11s} {'theta_stiff':>11s} {'Δθ_raw':>7s} "
        f"{'Δθ_corr':>8s} {'Φ_lane':>7s} {'n_s/n_st':>10s}")
    lines.append("-" * 90)
    for k in range(int(results["n_cells"])):
        s = float(results["G_s_soft"][k])
        raw = float(results['d_theta'][k])
        corr = raw - base.get(s, 0.0)
        lines.append(
            f"{results['G_s_soft'][k]:10.3f} "
            f"{results['G_s_stiff'][k]:10.3f} "
            f"{results['Ca_ratio'][k]:9.2f} "
            f"{results['theta_soft'][k]:+11.2f} "
            f"{results['theta_stiff'][k]:+11.2f} "
            f"{raw:7.2f} "
            f"{corr:+8.2f} "
            f"{results['lane_order'][k]:+7.2f} "
            f"{int(results['n_soft_placed'][k])}/"
            f"{int(results['n_stiff_placed'][k])}")
    lines.append("=" * 90)
    lines.append("")
    lines.append("Δθ_corr = Δθ_raw  −  Δθ_raw at the diagonal cell of the same row.")
    lines.append("Diagonal baselines:")
    for s, b in sorted(base.items()):
        lines.append(f"  G_s = {s:.3f}  →  Δθ_diag = {b:.2f}°")
    out_path.write_text("\n".join(lines) + "\n")


def main():
    npz_path = HERE / "sweep_results.npz"
    if not npz_path.is_file():
        print(f"Missing {npz_path}. Run sweep.py first.")
        return 2

    data = dict(np.load(npz_path, allow_pickle=False))
    write_text_summary(HERE / "fig_sweep_summary.txt", data)
    print(f"Wrote {HERE / 'fig_sweep_summary.txt'}")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; text summary only")
        return 0

    sweep_axis = np.asarray(data["sweep_values"])
    soft  = data["G_s_soft"]
    stiff = data["G_s_stiff"]
    dt    = data["d_theta"]
    lo    = data["lane_order"]
    ca    = data["Ca_ratio"]

    # ── Fig 1: Δθ heatmap on (G_s_soft, G_s_stiff) ─────────────────
    G_dtheta = _grid_from_pairs(soft, stiff, dt, sweep_axis)

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(G_dtheta, origin="lower", cmap="viridis",
                    aspect="equal")
    ax.set_xticks(range(sweep_axis.size))
    ax.set_xticklabels([f"{v:.2g}" for v in sweep_axis])
    ax.set_yticks(range(sweep_axis.size))
    ax.set_yticklabels([f"{v:.2g}" for v in sweep_axis])
    ax.set_xlabel("$G_s$  (stiff species)")
    ax.set_ylabel("$G_s$  (soft species)")
    ax.set_title("DLD displacement-angle separation $|Δθ|$  (deg)")
    # Overlay numerical values
    for i in range(sweep_axis.size):
        for j in range(sweep_axis.size):
            if np.isnan(G_dtheta[i, j]):
                continue
            ax.text(j, i, f"{G_dtheta[i, j]:.1f}",
                     ha="center", va="center",
                     color="white" if G_dtheta[i, j] >
                     np.nanmax(G_dtheta) / 2 else "black",
                     fontsize=9)
    fig.colorbar(im, ax=ax, label="|Δθ| (deg)")
    fig.tight_layout()
    fig.savefig(HERE / "fig_sweep_dtheta.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {HERE / 'fig_sweep_dtheta.png'}")

    # ── Fig 2: lane-order heatmap ──────────────────────────────────
    G_lo = _grid_from_pairs(soft, stiff, lo, sweep_axis)

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(G_lo, origin="lower", cmap="plasma",
                    aspect="equal", vmin=-1, vmax=1)
    ax.set_xticks(range(sweep_axis.size))
    ax.set_xticklabels([f"{v:.2g}" for v in sweep_axis])
    ax.set_yticks(range(sweep_axis.size))
    ax.set_yticklabels([f"{v:.2g}" for v in sweep_axis])
    ax.set_xlabel("$G_s$  (stiff species)")
    ax.set_ylabel("$G_s$  (soft species)")
    ax.set_title("Lane-order parameter  $\\Phi_{\\rm lane}$")
    for i in range(sweep_axis.size):
        for j in range(sweep_axis.size):
            if np.isnan(G_lo[i, j]):
                continue
            ax.text(j, i, f"{G_lo[i, j]:+.2f}",
                     ha="center", va="center", fontsize=8,
                     color="white" if abs(G_lo[i, j]) > 0.4 else "black")
    fig.colorbar(im, ax=ax, label="$\\Phi_{\\rm lane}$")
    fig.tight_layout()
    fig.savefig(HERE / "fig_sweep_lane_order.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {HERE / 'fig_sweep_lane_order.png'}")

    # ── Fig 3: Δθ vs Ca ratio (collapse onto one axis) ─────────────
    fig, ax = plt.subplots(figsize=(6, 4))
    finite = ~np.isnan(dt)
    ax.scatter(ca[finite], dt[finite], s=60, alpha=0.7,
                c=soft[finite], cmap="viridis")
    ax.set_xscale("log")
    ax.set_xlabel("Ca ratio $G_{s,\\,\\rm stiff} / G_{s,\\,\\rm soft}$")
    ax.set_ylabel("$|Δθ|$  (deg)")
    ax.set_title("Sorting power vs stiffness contrast")
    ax.axhline(0, color="gray", lw=0.5)
    ax.axhline(2.0, color="red", ls="--", lw=0.8,
                label="sorting threshold (2°)")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(HERE / "fig_sweep_dtheta_vs_Ca.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {HERE / 'fig_sweep_dtheta_vs_Ca.png'}")

    # ── Fig 4: diagonal-corrected Δθ heatmap ──────────────────────
    # Subtract the same-row diagonal value from each cell. The diagonal
    # cells (G_s_soft = G_s_stiff) should give Δθ ≈ 0 by construction;
    # any non-zero diagonal is a placement-induced systematic bias
    # that we subtract row-by-row to isolate the deformability signal.
    baseline = diagonal_baseline(data)
    dt_corr  = np.array([dt[k] - baseline.get(float(soft[k]), 0.0)
                          for k in range(len(dt))])
    G_corr   = _grid_from_pairs(soft, stiff, dt_corr, sweep_axis)

    fig, ax = plt.subplots(figsize=(6, 5))
    vmax = max(0.5, np.nanmax(np.abs(G_corr)))
    im = ax.imshow(G_corr, origin="lower", cmap="RdBu_r",
                    aspect="equal", vmin=-vmax, vmax=vmax)
    ax.set_xticks(range(sweep_axis.size))
    ax.set_xticklabels([f"{v:.2g}" for v in sweep_axis])
    ax.set_yticks(range(sweep_axis.size))
    ax.set_yticklabels([f"{v:.2g}" for v in sweep_axis])
    ax.set_xlabel("$G_s$  (stiff species)")
    ax.set_ylabel("$G_s$  (soft species)")
    ax.set_title("$Δθ$ minus same-row diagonal baseline  (deg)")
    for i in range(sweep_axis.size):
        for j in range(sweep_axis.size):
            if np.isnan(G_corr[i, j]):
                continue
            ax.text(j, i, f"{G_corr[i, j]:+.2f}",
                     ha="center", va="center", fontsize=9,
                     color="white" if abs(G_corr[i, j]) > 0.5 * vmax
                                       else "black")
    fig.colorbar(im, ax=ax, label="Δθ_corrected (deg)")
    fig.tight_layout()
    fig.savefig(HERE / "fig_sweep_dtheta_corrected.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {HERE / 'fig_sweep_dtheta_corrected.png'}")

    # ── Fig 5: corrected Δθ vs Ca ratio (THE HEADLINE) ────────────
    fig, ax = plt.subplots(figsize=(6, 4))
    finite_corr = ~np.isnan(dt_corr) & (np.abs(ca - 1.0) > 1e-9)  # drop diagonals
    sc = ax.scatter(ca[finite_corr], dt_corr[finite_corr],
                     s=70, alpha=0.85, edgecolors="black", lw=0.5,
                     c=soft[finite_corr], cmap="viridis")
    ax.set_xscale("log")
    ax.set_xlabel("Ca ratio $G_{s,\\,\\rm stiff} / G_{s,\\,\\rm soft}$")
    ax.set_ylabel("$Δθ_{\\rm corrected}$  (deg)")
    ax.set_title("Deformability-driven sorting (diagonal baseline removed)")
    ax.axhline(0, color="gray", lw=0.5)
    ax.grid(alpha=0.3)
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("$G_{s,\\,\\rm soft}$")
    fig.tight_layout()
    fig.savefig(HERE / "fig_sweep_dtheta_corrected_vs_Ca.png", dpi=160)
    plt.close(fig)
    print(f"Wrote {HERE / 'fig_sweep_dtheta_corrected_vs_Ca.png'}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
