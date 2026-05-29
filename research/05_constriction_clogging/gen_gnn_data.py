#!/usr/bin/env python3
"""Generate a REAL labeled graph dataset for the Flavor-B GNN.

Runs the constriction-clogging simulation (stable periodic + body-force path)
across a few apertures, extracts a contact graph from every output frame, and
labels each graph by its run's clog outcome. Graphs from all stable cells are
pooled into a single dataset that ``graph_gnn_demo.py`` (or any trainer) can
consume via ``pysoftflow.ml.load_graph_dataset``.

Diverged cells (NaN in the trajectory) are detected and skipped so one unstable
aperture cannot poison the dataset.

Usage:
    python research/05_constriction_clogging/gen_gnn_data.py
"""

from __future__ import annotations

import io
import contextlib
import pathlib
import sys

import numpy as np

script_dir = pathlib.Path(__file__).resolve().parent
project_dir = script_dir.parents[1]
sys.path.insert(0, str(project_dir / "python"))
for _d in ("build", "build_phase2", "build_phase1"):
    _c = project_dir / _d / "python"
    if _c.is_dir():
        sys.path.insert(0, str(_c))

import run as run_mod  # noqa: E402
from analyse import NX, NY, R_CAPSULE, analyse_trajectory  # noqa: E402
from pysoftflow.ml import build_graph_dataset, load_graph_dataset, save_graph_dataset  # noqa: E402

# Span narrow (clog) to wide (pass) with EXTRA sampling near the clog->pass
# transition (~D/d=3), using the DENSE deterministic "fill" placement (real
# contacts/arches). This single-seed, denser-aperture dataset gave the best
# honest generalization (LORO 0.871).
#
# NOTE (negative result): adding a jittered second seed per aperture for
# run-diversity was tried and *reduced* LORO to 0.70. Near the sharp clog/pass
# boundary the outcome is micro-config-sensitive, so two jittered seeds at the
# same aperture often land on OPPOSITE sides of the boundary -> label noise that
# hurts more than the near-twin helps. So we keep the deterministic fill.
APERTURES = [1.5, 2.0, 2.5, 2.75, 3.0, 3.25, 3.5, 4.0, 4.5, 5.0]
SEEDS = [0x05C106]
PLACEMENT = "fill"
FILL_JITTER = 0.0     # deterministic dense fill (jitter degraded generalization)
CA = 0.01
STEPS = 10000
FRAME_STRIDE = 2
NECK_X = 0.5 * (run_mod.NECK_X[0] + run_mod.NECK_X[1])
OUT = script_dir / "gnn_data"


def _trajectory_ok(csv: pathlib.Path) -> bool:
    """True if the trajectory exists and contains no NaN/Inf (no divergence)."""
    if not csv.exists():
        return False
    t = np.genfromtxt(csv, delimiter=",", names=True)
    return bool(np.all(np.isfinite(t["cx"])) and np.all(np.isfinite(t["cy"])))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pooled = []
    n_runs = len(APERTURES) * len(SEEDS)
    print(f"Generating GNN data: {len(APERTURES)} apertures x {len(SEEDS)} seeds "
          f"= {n_runs} runs ({PLACEMENT} placement, Ca={CA}, {STEPS} steps each)\n")

    for ap in APERTURES:
        for si, seed in enumerate(SEEDS):
            cell = OUT / f"ap{ap:.2f}_s{si}".replace(".", "p")
            cell.mkdir(parents=True, exist_ok=True)
            sim, meta = run_mod.build_sim(aperture=ap, ca=CA, output_dir=str(cell),
                                          rng_seed=seed, placement=PLACEMENT,
                                          jitter=FILL_JITTER)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                sim.thermo(interval=max(STEPS // 10, 200))
                sim.warmup(steps=500, ramp_steps=1500)
                sim.run(STEPS)
            csv = cell / "particle_data.csv"

            if not _trajectory_ok(csv):
                print(f"  D/d={ap:.2f} s{si}: DIVERGED or no output -> skipped")
                continue

            m = analyse_trajectory(str(csv), NECK_X)
            graphs = build_graph_dataset(
                str(csv), radius=R_CAPSULE, box=(NX, NY), periodic_x=True,
                contact_cutoff=1.0, frame_stride=FRAME_STRIDE, label=float(m.clogged),
            )
            save_graph_dataset(graphs, str(cell / "graphs.npz"))
            pooled.extend(graphs)
            print(f"  D/d={ap:.2f} s{si}: clogged={m.clogged}  "
                  f"escaped={m.escaped_fraction:.2f}  Z_neck={m.neck_contact_number:.2f}  "
                  f"-> {len(graphs)} graphs")

    if not pooled:
        print("\nNo stable cells produced data.")
        return

    save_graph_dataset(pooled, str(OUT / "gnn_dataset.npz"))
    n_clog = sum(1 for g in pooled if g.graph_label == 1.0)
    print(f"\nPooled REAL dataset: {len(pooled)} graphs "
          f"({n_clog} clog / {len(pooled) - n_clog} no-clog)")
    print(f"  -> {OUT / 'gnn_dataset.npz'}")
    # Sanity: reload + report contact statistics.
    re = load_graph_dataset(str(OUT / "gnn_dataset.npz"))
    deg = np.array([g.mean_degree for g in re])
    print(f"  reload OK: {len(re)} graphs, mean degree range "
          f"[{deg.min():.2f}, {deg.max():.2f}]")
    print("\nTrain on it with: from pysoftflow.ml import load_graph_dataset; "
          "from pysoftflow.ml.gnn import train_clog_gnn")


if __name__ == "__main__":
    main()
