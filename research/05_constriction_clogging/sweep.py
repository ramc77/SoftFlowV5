"""
Project 5 sweep — (D/d, Ca) clog phase diagram.
================================================
Runs the constriction-clogging simulation across a grid of aperture ratios
D/d and capillary numbers Ca, recording the clog flag and escapee metrics for
each cell. The output is the clog/no-clog phase diagram to compare against the
Bielinski 2021 state diagram (rigid clog threshold D/d ~ 3; soft escape for
Ca > ~0.005) and the 2-D Hong 2017 / Marin 2018 experiments.

Each cell also dumps a GNN graph dataset (one graph per frame, labelled by the
cell's clog outcome) so the full sweep doubles as a training set for the
Flavor-B graph neural network.

Usage:
    python research/05_constriction_clogging/sweep.py            # full grid
    python research/05_constriction_clogging/sweep.py --smoke    # 5k steps/cell
"""

from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np

script_dir = pathlib.Path(__file__).resolve().parent
project_dir = script_dir.parents[1]
sys.path.insert(0, str(project_dir / "python"))

import run as run_mod  # noqa: E402  (project-5 run.py)
from analyse import analyse_trajectory  # noqa: E402

# Sweep axes (Bielinski 2021 spanned D/d in [1.25, 3.5], Ca in [1e-3, 1e-1]).
APERTURES = (1.25, 1.5, 2.0, 2.5, 3.0, 3.5)
CA_VALUES = (1e-3, 1e-2, 5e-2, 1e-1)

NECK_X = 0.5 * (run_mod.NECK_X[0] + run_mod.NECK_X[1])


def run_cell(aperture: float, ca: float, n_steps: int, out_root: pathlib.Path) -> dict:
    """Run one (D/d, Ca) cell and return its result record."""
    cell_dir = out_root / f"ap{aperture:.2f}_ca{ca:.0e}".replace(".", "p")
    cell_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    sim, meta = run_mod.build_sim(aperture=aperture, ca=ca, output_dir=str(cell_dir))
    sim.thermo(interval=max(n_steps // 20, 100))
    sim.warmup(steps=500, ramp_steps=1500)
    sim.run(n_steps)
    walltime = time.time() - t0

    csv_path = cell_dir / "particle_data.csv"
    metrics = analyse_trajectory(
        str(csv_path), NECK_X, save_graphs=str(cell_dir / "graphs.npz")
    )

    record = {
        "aperture": aperture,
        "ca": ca,
        "neck_width": meta["neck_width"],
        "g_s": meta["g_s"],
        "clogged": int(metrics.clogged),
        "escaped_fraction": metrics.escaped_fraction,
        "evacuation_time": metrics.evacuation_time,
        "neck_contact_number": metrics.neck_contact_number,
        "walltime_s": walltime,
        "_n_steps": n_steps,
    }
    (cell_dir / "cell_result.json").write_text(json.dumps(record, indent=2))
    return record


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--smoke", action="store_true", help="5000 steps per cell")
    p.add_argument("--steps", type=int, default=None)
    p.add_argument("--out", type=str, default=str(script_dir / "sweep_out"))
    args = p.parse_args(argv if argv is not None else sys.argv[1:])

    n_steps = args.steps if args.steps is not None else (5_000 if args.smoke else 60_000)
    out_root = pathlib.Path(args.out)
    out_root.mkdir(parents=True, exist_ok=True)

    records = []
    n_total = len(APERTURES) * len(CA_VALUES)
    for k, ap in enumerate(APERTURES):
        for j, ca in enumerate(CA_VALUES):
            idx = k * len(CA_VALUES) + j + 1
            print(f"[{idx}/{n_total}] D/d={ap:.2f}, Ca={ca:.0e} ...", flush=True)
            records.append(run_cell(ap, ca, n_steps, out_root))

    # Aggregate to a flat npz (mirrors projects 03/04).
    keys = ("aperture", "ca", "clogged", "escaped_fraction",
            "evacuation_time", "neck_contact_number", "g_s", "walltime_s")
    agg = {key: np.array([r[key] for r in records], dtype=float) for key in keys}
    agg["aperture_values"] = np.array(APERTURES, dtype=float)
    agg["ca_values"] = np.array(CA_VALUES, dtype=float)
    agg["n_cells"] = np.array(n_total)
    np.savez(out_root / "sweep_results.npz", **agg)
    print(f"\nWrote {out_root / 'sweep_results.npz'}  ({n_total} cells)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
