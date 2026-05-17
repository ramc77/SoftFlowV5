"""Plot Project-3 results from history.npz.

Three publication-figure candidates:

  fig1_trajectories.png   — overlaid particle trajectories, coloured by
                              type, against the pillar-array sketch.
  fig2_dld_angle.png      — distribution of per-particle DLD angle
                              θ_DLD, soft vs stiff.
  fig3_summary_bars.png   — mean θ_DLD, lane order, and the Z_ij
                              contact matrix as a heatmap.
"""

import os
import pathlib
import sys

import numpy as np


def main():
    here = pathlib.Path(__file__).parent
    h    = np.load(here / "history.npz")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; skipping plots")
        return 0

    traj_xy = h["traj_xy"]            # (n_snap, n_part, 2)
    types_  = h["types"]
    dx      = h["dx_per_p"]
    dy      = h["dy_per_p"]

    soft_mask  = types_ == 0
    stiff_mask = types_ == 1

    # ── Fig 1: trajectories ────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 4))
    for i in range(traj_xy.shape[1]):
        col = "#1f77b4" if soft_mask[i] else "#d62728"
        # Plot only paths that did not wrap (continuous in x).
        x_path = traj_xy[:, i, 0]
        y_path = traj_xy[:, i, 1]
        # Mask jumps from wrap.
        dxstep = np.diff(x_path)
        wraps  = np.where(dxstep < -50)[0]
        if len(wraps) == 0:
            ax.plot(x_path, y_path, color=col, lw=0.7, alpha=0.6)
        else:
            # Plot in segments split at each wrap.
            starts = np.concatenate(([0], wraps + 1))
            ends   = np.concatenate((wraps + 1, [len(x_path)]))
            for s, e in zip(starts, ends):
                ax.plot(x_path[s:e], y_path[s:e],
                         color=col, lw=0.7, alpha=0.6)
    ax.set_xlabel("x (lattice)")
    ax.set_ylabel("y (lattice)")
    ax.set_title("Particle trajectories (blue = soft, red = stiff)")
    fig.tight_layout()
    fig.savefig(here / "fig1_trajectories.png", dpi=160)
    plt.close(fig)

    # ── Fig 2: DLD-angle distribution ──────────────────────────────
    fig, ax = plt.subplots(figsize=(6, 4))
    theta_p = np.degrees(np.arctan2(dy, dx))
    ax.hist([theta_p[soft_mask], theta_p[stiff_mask]],
             bins=30, label=["soft", "stiff"],
             color=["#1f77b4", "#d62728"], alpha=0.7,
             histtype="stepfilled", edgecolor="black", lw=0.5)
    ax.axvline(float(h["theta_soft"]),  ls="--", color="#1f77b4",
                label=f"<θ>_soft  = {float(h['theta_soft']):+.2f}°")
    ax.axvline(float(h["theta_stiff"]), ls="--", color="#d62728",
                label=f"<θ>_stiff = {float(h['theta_stiff']):+.2f}°")
    ax.set_xlabel("DLD displacement angle  θ (deg)")
    ax.set_ylabel("count")
    ax.set_title("Per-particle DLD angle distribution")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(here / "fig2_dld_angle.png", dpi=160)
    plt.close(fig)

    # ── Fig 3: summary + contact matrix ────────────────────────────
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.5))

    # mean θ_DLD bars
    ax[0].bar(["soft", "stiff"],
                [float(h["theta_soft"]), float(h["theta_stiff"])],
                color=["#1f77b4", "#d62728"])
    ax[0].axhline(0, color="black", lw=0.5)
    ax[0].set_ylabel("<θ_DLD>  (deg)")
    ax[0].set_title("DLD angle per type")

    # lane order parameter
    ax[1].bar(["lane order"],
                [float(h["lane_order"])], color="#2ca02c")
    ax[1].set_ylim(-1.0, 1.0)
    ax[1].axhline(0, color="black", lw=0.5)
    ax[1].set_ylabel("Φ_lane")
    ax[1].set_title("Lane order (x-axis)")

    # contact matrix heatmap
    Z = h["Z_matrix"]
    im = ax[2].imshow(Z, cmap="viridis", vmin=0)
    ax[2].set_xticks([0, 1]); ax[2].set_xticklabels(["soft", "stiff"])
    ax[2].set_yticks([0, 1]); ax[2].set_yticklabels(["soft", "stiff"])
    for i in range(2):
        for j in range(2):
            ax[2].text(j, i, f"{Z[i, j]:.2f}",
                        ha="center", va="center",
                        color="white" if Z[i, j] > Z.max() / 2 else "black")
    ax[2].set_title("Z_ij per-type contact matrix")
    fig.colorbar(im, ax=ax[2], shrink=0.7)

    fig.tight_layout()
    fig.savefig(here / "fig3_summary_bars.png", dpi=160)
    plt.close(fig)

    print(f"Wrote 3 figures to {here}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
