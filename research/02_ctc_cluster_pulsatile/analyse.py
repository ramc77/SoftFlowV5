"""Plot Project-2 results from history.npz.

Three publication-figure candidates:

  fig1_cluster_size.png  — largest cluster vs time, steady vs pulsatile
  fig2_bonds_vs_force.png — bond count phase-portrait (bonds vs body_fx)
  fig3_summary_bars.png   — peak/mean cluster size + bonds (the headline)
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

    steady   = {k.replace("steady_", ""):  h[k]
                  for k in h.files if k.startswith("steady_")}
    pulsatile = {k.replace("pulsatile_", ""): h[k]
                   for k in h.files if k.startswith("pulsatile_")}

    # ── Fig 1: cluster size vs time ────────────────────────────────
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(steady["time"],   steady["largest_cluster"],
             label="steady",    color="#1f77b4", lw=1.5)
    ax.plot(pulsatile["time"], pulsatile["largest_cluster"],
             label="pulsatile", color="#d62728", lw=1.5)
    ax.set_xlabel("time (lattice)")
    ax.set_ylabel("largest cluster size (cells)")
    ax.set_title("CTC cluster size vs time")
    ax.legend()
    fig.tight_layout()
    fig.savefig(here / "fig1_cluster_size.png", dpi=160)
    plt.close(fig)

    # ── Fig 2: bonds vs body force (phase portrait) ────────────────
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(steady["body_fx"], steady["n_bonds"],
                s=10, alpha=0.3, color="#1f77b4", label="steady")
    ax.scatter(pulsatile["body_fx"], pulsatile["n_bonds"],
                s=10, alpha=0.4, color="#d62728", label="pulsatile")
    ax.set_xlabel("body force F_x  (lattice)")
    ax.set_ylabel("active bond count")
    ax.set_title("Bond population vs instantaneous driving force")
    ax.legend()
    fig.tight_layout()
    fig.savefig(here / "fig2_bonds_vs_force.png", dpi=160)
    plt.close(fig)

    # ── Fig 3: summary bars ───────────────────────────────────────
    fig, ax = plt.subplots(1, 2, figsize=(9, 4))

    labels = ["steady", "pulsatile"]
    peaks  = [steady["largest_cluster"].max(),
                pulsatile["largest_cluster"].max()]
    means  = [steady["largest_cluster"].mean(),
                pulsatile["largest_cluster"].mean()]
    bonds  = [steady["n_bonds"].mean(),
                pulsatile["n_bonds"].mean()]

    ax[0].bar(labels, peaks, color=["#1f77b4", "#d62728"])
    ax[0].set_ylabel("peak cluster size")
    ax[0].set_title("Peak cluster size")
    for i, v in enumerate(peaks):
        ax[0].text(i, v, f"{int(v)}", ha="center", va="bottom")

    ax[1].bar(labels, bonds, color=["#1f77b4", "#d62728"])
    ax[1].set_ylabel("mean active bonds")
    ax[1].set_title("Mean bond count")
    for i, v in enumerate(bonds):
        ax[1].text(i, v, f"{v:.1f}", ha="center", va="bottom")

    fig.tight_layout()
    fig.savefig(here / "fig3_summary_bars.png", dpi=160)
    plt.close(fig)

    print(f"Wrote 3 figures to {here}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
