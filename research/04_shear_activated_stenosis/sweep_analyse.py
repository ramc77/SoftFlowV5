"""
Project 4 — phase-diagram producer
==================================
Reads ``sweep_results.npz`` and writes:

    fig_eta_deposit.png         [HEADLINE]  η_deposit × (severity, γ̇_th)
    fig_eta_offtarget.png       off-target loss
    fig_release_fraction.png    carrier release activity
    fig_t50_deposit.png         log10(t_50_deposit)
    fig_sweet_spot.png          η_deposit vs γ̇_th per severity (story plot)
    fig_sweep_summary.txt       tabular summary

The sweet-spot ridge is the publishable design rule:
γ̇_th^star(σ) ≈ throat shear rate at severity σ.
"""

from __future__ import annotations

import pathlib
import sys
import numpy as np
import matplotlib.pyplot as plt


script_dir = pathlib.Path(__file__).resolve().parent
RESULTS_FILE = script_dir / "sweep_results.npz"


def load_results():
    if not RESULTS_FILE.exists():
        sys.exit(f"sweep_results.npz not found at {RESULTS_FILE}.\n"
                  f"Run sweep.py first.")
    data = np.load(RESULTS_FILE, allow_pickle=True)
    cols = {k: data[k] for k in data.files
             if k not in ("severity_values", "gamma_th_values",
                          "n_steps", "n_cells")}
    sev_axis = np.asarray(data["severity_values"], dtype=int)
    th_axis  = np.asarray(data["gamma_th_values"], dtype=float)
    return cols, sev_axis, th_axis


def reshape_grid(values, sev_col, th_col, sev_axis, th_axis,
                  default=np.nan):
    n_s, n_t = len(sev_axis), len(th_axis)
    grid = np.full((n_s, n_t), default, dtype=float)
    sev_col = np.asarray(sev_col, dtype=int)
    th_col  = np.asarray(th_col, dtype=float)
    for v, s, t in zip(values, sev_col, th_col):
        i = int(np.argmin(np.abs(sev_axis - s)))
        j = int(np.argmin(np.abs(th_axis - t)))
        grid[i, j] = v
    return grid


def _heatmap(ax, grid, sev_axis, th_axis, *,
              cmap="RdPu", title="", cbar_label="",
              fmt="{:.3f}", vmin=None, vmax=None):
    im = ax.imshow(grid, origin="lower", aspect="auto",
                    cmap=cmap, vmin=vmin, vmax=vmax)
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            v = grid[i, j]
            label = "—" if np.isnan(v) else fmt.format(v)
            ax.text(j, i, label, ha="center", va="center",
                    fontsize=9, color="black")
    ax.set_xticks(range(len(th_axis)))
    ax.set_yticks(range(len(sev_axis)))
    ax.set_xticklabels([f"{t:.0e}" for t in th_axis], rotation=30, ha="right")
    ax.set_yticklabels([f"{s}%" for s in sev_axis])
    ax.set_xlabel(r"Shear activation threshold $\dot\gamma_\mathrm{th}$ (LU)")
    ax.set_ylabel(r"Stenosis severity $\sigma$")
    ax.set_title(title)
    cb = plt.colorbar(im, ax=ax)
    cb.set_label(cbar_label)


def fig_eta_deposit(cols, sev_axis, th_axis, out_path):
    """HEADLINE: deposition efficiency phase diagram."""
    grid = reshape_grid(cols["eta_deposit"],
                         cols["severity"], cols["gamma_th"], sev_axis, th_axis)
    fig, ax = plt.subplots(figsize=(7, 5))
    vmax = max(0.05, min(1.0, float(np.nanmax(grid))))
    _heatmap(ax, grid, sev_axis, th_axis,
              cmap="RdPu",
              title=r"Deposition efficiency $\eta_\mathrm{deposit}$",
              cbar_label=r"$\eta_\mathrm{deposit}$",
              fmt="{:.3f}", vmin=0.0, vmax=vmax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def fig_eta_offtarget(cols, sev_axis, th_axis, out_path):
    grid = reshape_grid(cols["eta_offtarget"],
                         cols["severity"], cols["gamma_th"], sev_axis, th_axis)
    fig, ax = plt.subplots(figsize=(7, 5))
    vmax = max(0.05, min(1.0, float(np.nanmax(grid))))
    _heatmap(ax, grid, sev_axis, th_axis,
              cmap="Blues",
              title=r"Off-target fraction $\eta_\mathrm{offtarget}$",
              cbar_label=r"$\eta_\mathrm{offtarget}$",
              fmt="{:.3f}", vmin=0.0, vmax=vmax)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def fig_release_fraction(cols, sev_axis, th_axis, out_path):
    grid = reshape_grid(cols["release_fraction"],
                         cols["severity"], cols["gamma_th"], sev_axis, th_axis)
    fig, ax = plt.subplots(figsize=(7, 5))
    _heatmap(ax, grid, sev_axis, th_axis,
              cmap="YlOrBr",
              title=r"Carrier release fraction $M_\mathrm{rel}/M_\mathrm{loaded}$",
              cbar_label="released / loaded",
              fmt="{:.3f}", vmin=0.0, vmax=1.0)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def fig_t50(cols, sev_axis, th_axis, out_path):
    grid = reshape_grid(cols["t_50_deposit"],
                         cols["severity"], cols["gamma_th"], sev_axis, th_axis)
    with np.errstate(invalid="ignore", divide="ignore"):
        logt = np.log10(grid)
    fig, ax = plt.subplots(figsize=(7, 5))
    _heatmap(ax, logt, sev_axis, th_axis,
              cmap="PuBu",
              title=r"$\log_{10}$ time-to-50%-deposit",
              cbar_label=r"$\log_{10} t_{50}$",
              fmt="{:.2f}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def fig_sweet_spot(cols, sev_axis, th_axis, out_path):
    """η_deposit vs γ̇_th, one curve per severity. The sweet-spot ridge."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    sev_col = np.asarray(cols["severity"], dtype=int)
    th_col  = np.asarray(cols["gamma_th"], dtype=float)
    eta_col = np.asarray(cols["eta_deposit"])
    for sev in sev_axis:
        mask = sev_col == int(sev)
        if not mask.any():
            continue
        ths = th_col[mask]
        etas = eta_col[mask]
        order = np.argsort(ths)
        ax.plot(ths[order], etas[order], "o-", linewidth=1.7, markersize=7,
                label=fr"$\sigma = {sev}\%$")
    ax.set_xscale("log")
    ax.set_xlabel(r"Shear threshold $\dot\gamma_\mathrm{th}$ (LU)")
    ax.set_ylabel(r"$\eta_\mathrm{deposit}$")
    ax.set_title("Deposition efficiency vs threshold, by stenosis severity")
    ax.legend(loc="best", fontsize=9, frameon=True)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=-0.02)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def write_summary(cols, sev_axis, th_axis, out_path):
    lines = []
    lines.append("Project 4 — Shear-activated stenosis sweep summary")
    lines.append("=" * 100)
    lines.append("  sev   γ̇_th       η_deposit  η_offtarget  rel_frac   t_50_deposit   walltime")
    lines.append("-" * 100)
    pairs = sorted(zip(cols["severity"], cols["gamma_th"],
                        cols["eta_deposit"], cols["eta_offtarget"],
                        cols["release_fraction"],
                        cols["t_50_deposit"], cols["walltime_s"]),
                    key=lambda r: (r[0], r[1]))
    for (sev, th, ed, eo, rf, t50, wt) in pairs:
        ed_s  = f"{ed:.3f}"  if not np.isnan(ed)  else "  —  "
        eo_s  = f"{eo:.3f}"  if not np.isnan(eo)  else "  —  "
        t50_s = f"{t50:7.0f}" if not np.isnan(t50) else "    —  "
        lines.append(f"  {int(sev):>3d}%  {th:.1e}    {ed_s}      {eo_s}     "
                      f"{rf:.3f}      {t50_s}    {wt:6.1f}s")
    lines.append("=" * 100)
    lines.append("")
    # Phase classification
    lines.append("Phase classification by η_deposit:")
    for (sev, th, ed, *_) in pairs:
        if np.isnan(ed):
            tag = "NO RELEASE (γ̇_th too high)"
        elif ed < 0.05:
            tag = "FAILED DELIVERY"
        elif ed < 0.2:
            tag = "WEAK DELIVERY"
        elif ed < 0.5:
            tag = "MODERATE DELIVERY"
        else:
            tag = "STRONG DELIVERY"
        lines.append(f"  σ={int(sev):>3d}%  γ̇_th={th:.1e}  η={ed if not np.isnan(ed) else 0:.3f}  {tag}")
    out_path.write_text("\n".join(lines))
    print(f"Wrote {out_path}")


def main(argv=None) -> int:
    cols, sev_axis, th_axis = load_results()
    fig_eta_deposit(cols, sev_axis, th_axis,
                     script_dir / "fig_eta_deposit.png")
    fig_eta_offtarget(cols, sev_axis, th_axis,
                       script_dir / "fig_eta_offtarget.png")
    fig_release_fraction(cols, sev_axis, th_axis,
                          script_dir / "fig_release_fraction.png")
    fig_t50(cols, sev_axis, th_axis,
             script_dir / "fig_t50_deposit.png")
    fig_sweet_spot(cols, sev_axis, th_axis,
                    script_dir / "fig_sweet_spot.png")
    write_summary(cols, sev_axis, th_axis,
                   script_dir / "fig_sweep_summary.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
