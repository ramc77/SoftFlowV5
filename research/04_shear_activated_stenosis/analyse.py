"""
Project 4 — single-cell deep-dive analyser
==========================================
Reads ``history.npz`` and emits a 4-panel figure:

  A — cumulative drug released by carriers
  B — cumulative deposition at the stenosis (TARGET)
  C — cumulative loss to off-target walls
  D — instantaneous deposition efficiency η_deposit(t)

Usage
-----
    python analyse.py                        # ./history.npz
    python analyse.py sweep_out/sev50_th2/   # any per-cell directory
"""

from __future__ import annotations

import pathlib
import sys
import numpy as np
import matplotlib.pyplot as plt


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    cell_dir = (pathlib.Path(argv[0]).resolve() if argv
                 else pathlib.Path(__file__).resolve().parent)
    f = cell_dir / "history.npz"
    if not f.exists():
        sys.exit(f"history.npz not found in {cell_dir}.")
    data = np.load(f, allow_pickle=True)
    h = data["history"]
    severity = int(data["severity"]) if "severity" in data.files else 0
    gamma_th = float(data["gamma_th"]) if "gamma_th" in data.files else 0.0
    t = h["time"]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)

    ax = axes[0, 0]
    ax.plot(t, h["total_released"], color="#2c3e50", lw=1.6)
    ax.set_ylabel("released")
    ax.set_title("A — Cumulative payload released by carriers")
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.plot(t, h["target_absorbed"], color="#c0392b", lw=1.6)
    ax.set_ylabel("target uptake")
    ax.set_title("B — Cumulative drug absorbed at stenosis")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    ax.plot(t, h["off_target_absorbed"], color="#7f8c8d", lw=1.6)
    ax.set_xlabel("step")
    ax.set_ylabel("off-target uptake")
    ax.set_title("C — Cumulative drug lost off-target")
    ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    rel = h["total_released"]
    upt = h["target_absorbed"]
    with np.errstate(invalid="ignore", divide="ignore"):
        eta = np.where(rel > 1e-12, upt / rel, 0.0)
    ax.plot(t, eta, color="#16a085", lw=1.6)
    last_q = int(0.75 * len(eta))
    eta_final = float(np.nanmean(eta[last_q:])) if last_q < len(eta) else 0.0
    ax.axhline(eta_final, color="red", linestyle="--", alpha=0.6,
                label=fr"$\eta_\mathrm{{deposit, final}} = {eta_final:.3f}$")
    ax.set_xlabel("step")
    ax.set_ylabel(r"$\eta_\mathrm{deposit}(t)$")
    ax.set_title(r"D — Deposition efficiency $\eta_\mathrm{deposit}$")
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=-0.02)

    fig.suptitle(fr"Cell: $\sigma = {severity}\%$,  "
                  fr"$\dot\gamma_\mathrm{{th}} = {gamma_th:.1e}$",
                  fontsize=11, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = cell_dir / "fig_cell_anatomy.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
