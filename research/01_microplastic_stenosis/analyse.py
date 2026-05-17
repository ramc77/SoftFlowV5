"""Plot Project-1 results from history.npz.

Produces three publication-figure candidates:

  fig1_dose_map.png       — 2D heatmap of leached chemical concentration
                              with stenosis outlined.
  fig2_per_particle.png   — M_p remaining vs particle type (boxplot or
                              swarm) showing size-dependent leaching.
  fig3_wall_vs_bulk.png   — Bar chart: dose at stenotic walls vs bulk
                              (the headline number for the paper).

Skipped silently if matplotlib isn't installed.
"""

import os
import pathlib
import sys

import numpy as np


def main():
    here = pathlib.Path(__file__).parent
    h = np.load(here / "history.npz")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle
    except ImportError:
        print("matplotlib not installed; skipping plots")
        print(f"  scalar field shape: {h['scalar_field'].shape}")
        print(f"  dose_bottom_at_stenosis: {h['dose_bottom_at_stenosis']}")
        print(f"  dose_top_at_stenosis:    {h['dose_top_at_stenosis']}")
        print(f"  dose_bulk:               {h['dose_bulk']}")
        return 0

    C = h["scalar_field"]
    ny, nx = C.shape
    sten_x = int(h["stenosis_x"])

    # ── Fig 1: spatial dose map ─────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 3.5))
    im = ax.imshow(C, origin="lower", aspect="auto",
                    extent=[0, nx, 0, ny],
                    cmap="inferno",
                    vmax=np.percentile(C, 99))
    ax.add_patch(Circle((sten_x, 8),       15, fill=False, color="white", lw=1.5))
    ax.add_patch(Circle((sten_x, ny - 8),  15, fill=False, color="white", lw=1.5))
    ax.set_xlabel("x (lattice)")
    ax.set_ylabel("y (lattice)")
    ax.set_title("Leached chemical concentration at end of run")
    fig.colorbar(im, ax=ax, label="C (lattice)")
    fig.tight_layout()
    fig.savefig(here / "fig1_dose_map.png", dpi=160)
    plt.close(fig)

    # ── Fig 2: per-particle M_p remaining by type ─────────────────
    labels  = ["RBC", "WBC", "MP small", "MP large"]
    M_p     = h["M_p_remaining"]
    type_   = h["type_label"]
    by_type = [M_p[type_ == t] for t in range(4)]

    fig, ax = plt.subplots(figsize=(6, 4))
    parts = ax.boxplot(by_type, labels=labels, showmeans=True, widths=0.4)
    ax.set_ylabel("Remaining particle mass M_p")
    ax.set_title("M_p reservoir per type after run")
    fig.tight_layout()
    fig.savefig(here / "fig2_per_particle.png", dpi=160)
    plt.close(fig)

    # ── Fig 3: wall vs bulk dose (the headline) ─────────────────────
    fig, ax = plt.subplots(figsize=(5, 4))
    db   = float(h["dose_bottom_at_stenosis"])
    dt   = float(h["dose_top_at_stenosis"])
    dbk  = float(h["dose_bulk"])
    bars = ["Bottom\nwall", "Top\nwall", "Bulk"]
    vals = [db, dt, dbk]
    cols = ["#d62728", "#d62728", "#1f77b4"]
    ax.bar(bars, vals, color=cols)
    ax.set_ylabel("Cumulative dose (lattice)")
    ax.set_title("Wall-localised leaching at stenosis")
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:.2g}", ha="center", va="bottom")
    fig.tight_layout()
    fig.savefig(here / "fig3_wall_vs_bulk.png", dpi=160)
    plt.close(fig)

    print(f"Wrote 3 figures to {here}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
