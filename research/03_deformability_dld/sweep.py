"""
Parameter sweep over (G_s_soft, G_s_stiff) for the DLD-deformability paper
==========================================================================
Runs the Project-3 simulation across a grid of stiffness pairs and
records the DLD displacement-angle separation Δθ, lane order Φ_lane,
and per-type contact matrix Z_ij for each cell.

Sweep grid
----------
Default: 5-value grid with the *upper triangle* (G_s_soft <
G_s_stiff), 10 cells total. The diagonal (G_s_soft = G_s_stiff) is
the no-contrast control where Δθ → 0 is expected. Skipping the
lower triangle saves time and avoids re-running mirror cells.

Output layout
-------------
research/03_deformability_dld/
├── sweep_results.npz          aggregated grid: theta_soft, theta_stiff,
│                                Δθ, lane_order, Z_matrix, n_placed
├── sweep_out/
│   ├── G0.030_S0.300/         per-cell directory (one for each pair)
│   │   ├── vtk_dld/...           VTK trajectory
│   │   ├── history.npz           per-cell time series
│   │   └── config/run_manifest.json  Phase-1 provenance
│   ├── G0.030_S0.600/
│   └── ...
└── sweep_analyse.py            heatmap producer (run me next)

Usage
-----
    python sweep.py --smoke         # 3-cell quick sanity (~3 min)
    python sweep.py                  # 10-cell upper-triangle (~30 min smoke / ~3 h prod)
    python sweep.py --full           # full 5×5 grid (~5 h prod)

Each cell writes its own run_manifest.json (git SHA, compiler flags,
RNG seed, fully resolved params) so every cell of the sweep is
independently reviewer-traceable.
"""

from __future__ import annotations

import argparse
import os
import pathlib
import random
import sys
import time
from dataclasses import dataclass

import numpy as np


# -- Path setup (find the build) --
script_dir   = pathlib.Path(__file__).resolve().parent
project_dir  = script_dir.parents[1]
sys.path.insert(0, str(project_dir / "python"))
for d in ("build", "build_phase1", "build_phase2"):
    cand = project_dir / d / "python"
    if cand.is_dir():
        sys.path.insert(0, str(cand))

from pysoftflow import SoftFlowSimulation                          # noqa: E402
from pysoftflow.analysis import SimulationSnapshot                 # noqa: E402
from pysoftflow.analysis.jamming import per_type_contact_stats     # noqa: E402
from pysoftflow.analysis.patterns import lane_order                # noqa: E402


# =================================================================
# Sweep grid + geometry (geometry identical to run.py)
# =================================================================

SWEEP_VALUES_FULL  = [0.03, 0.10, 0.20, 0.40, 0.80]   # 5-value grid; max 0.80 stable at MAX_LATTICE_FORCE=0.04
SWEEP_VALUES_SMOKE = [0.03, 0.10, 0.30]               # 3-value sub-grid (all known-safe)

NX, NY            = 400, 80
TAU               = 0.7
BODY_FX           = 3e-5
MAX_LATTICE_FORCE = 0.04

# DLD pillar array (must match run.py exactly so trajectory comparisons hold).
N_PILLAR_COLS  = 8
N_PILLAR_ROWS  = 4
PILLAR_R       = 4.0
PILLAR_DX      = 35.0
PILLAR_DY      = 14.0
ROW_SHIFT      = 3.5
ARRAY_X_START  = 50.0
ARRAY_Y_CENTRE = NY / 2.0

N_SOFT     = 20
N_STIFF    = 20
R_PART     = 3.0
OUT_EVERY  = 25

WARMUP_STEPS = 500     # shorter than the standalone example to keep sweep wall-time manageable
WARMUP_RAMP  = 1500


# =================================================================
# Result dataclass
# =================================================================


@dataclass
class CellResult:
    """One cell of the sweep — what we accumulate into sweep_results.npz."""
    G_s_soft:      float
    G_s_stiff:     float
    Ca_ratio:      float                # G_s_stiff / G_s_soft
    theta_soft:    float                # DLD angle, deg
    theta_stiff:   float
    d_theta:       float                # |Δθ|, deg (the headline)
    mean_dx_soft:  float
    mean_dy_soft:  float
    mean_dx_stiff: float
    mean_dy_stiff: float
    lane_order:    float
    Z00:           float                # soft–soft
    Z01:           float                # soft–stiff
    Z10:           float                # stiff–soft
    Z11:           float                # stiff–stiff
    n_soft_placed: int
    n_stiff_placed: int
    walltime_s:    float


# =================================================================
# One simulation cell
# =================================================================


def run_one_cell(*,
                  G_s_soft: float,
                  G_s_stiff: float,
                  output_dir: str,
                  n_steps: int,
                  rng_seed: int = 0x0AB10001) -> CellResult:
    """Build a DLD simulation for one (G_s_soft, G_s_stiff) pair, run
    it, return the aggregated CellResult."""
    t_start = time.time()

    sim = SoftFlowSimulation()
    sim.domain(nx=NX, ny=NY)
    sim.boundary(x="periodic", y="wall")
    sim.fluid(tau=TAU, max_lattice_force=MAX_LATTICE_FORCE,
               collision="regularized")
    sim.body_force(BODY_FX, 0.0)

    # DLD pillar array — identical to run.py.
    for col in range(N_PILLAR_COLS):
        for row in range(N_PILLAR_ROWS):
            cx = ARRAY_X_START + col * PILLAR_DX
            cy = (ARRAY_Y_CENTRE
                  + (row - (N_PILLAR_ROWS - 1) / 2.0) * PILLAR_DY
                  + col * ROW_SHIFT)
            if cy < 6 or cy > NY - 6:
                continue
            sim.obstacle("circle", center=(cx, cy), radius=PILLAR_R)

    sim.ibm(iterations=2)
    sim.lubrication(enabled=True)

    # Particle types — stiffness G_s is the swept parameter. All
    # other membrane params (C_skalak, k_bend, k_area, k_perimeter)
    # are held identical between the two types so the diagonal cell
    # (G_s_soft = G_s_stiff) is a true no-contrast control where
    # Δθ ≡ 0 by symmetry.
    K_BEND      = 0.010
    K_AREA      = 0.75
    K_PERIMETER = 0.075
    sim.particle_type("soft",
                       model="skalak",
                       G_s=G_s_soft,
                       C_skalak=10.0,
                       k_bend=K_BEND, k_area=K_AREA, k_perimeter=K_PERIMETER)
    sim.particle_type("stiff",
                       model="skalak",
                       G_s=G_s_stiff,
                       C_skalak=10.0,
                       k_bend=K_BEND, k_area=K_AREA, k_perimeter=K_PERIMETER)

    # Seeding bias fix (vs the original sweep). Two-pass RSA in the
    # same region biases the second species (first claims best
    # positions). Disjoint x-staggered inlets bias Δθ via different
    # x-entry into the periodic pillar array. Strictly alternating
    # singletons (stiff, soft, stiff, soft, …) still biases because
    # stiff always sees one fewer existing particle than soft at the
    # same step.
    #
    # Robust fix: **randomized-order singleton seeding** into a single
    # shared inlet. Build a list with N_STIFF "stiff" + N_SOFT "soft"
    # entries, shuffle with the sweep's RNG seed, and place one at a
    # time. Each species sees the same expected density at placement
    # time, and the species labels are statistically interchangeable.
    sim.region("inlet",
                x=(5, 47),
                y=(8, NY - 8))
    seeding_order = (["stiff"] * N_STIFF) + (["soft"] * N_SOFT)
    seeding_rng = random.Random(rng_seed)
    seeding_rng.shuffle(seeding_order)
    for k, species in enumerate(seeding_order):
        sim.generate(species,
                      count=1, region="inlet",
                      radius=R_PART, num_nodes=18,
                      method="random",
                      seed=42 + k,
                      min_gap=1.5)

    sim.output(format="vtk",
                directory=os.path.join(output_dir, "vtk_dld"),
                interval=OUT_EVERY)

    sim.initialize()

    # Trajectory logger — only the first and last snapshots are needed
    # for Δθ, but we save the full series in case the analysis script
    # wants to plot trajectories per cell.
    trajectory: list[dict] = []

    def _log(core, step):
        if step % OUT_EVERY != 0 and step != 0:
            return
        snap = SimulationSnapshot.from_simulation(core)
        trajectory.append({
            "step": int(step),
            "positions": snap.positions.copy(),
            "types":     snap.types.copy(),
        })

    sim.core.setStepCallback(_log)
    sim.warmup(steps=WARMUP_STEPS, ramp_steps=WARMUP_RAMP)
    sim.thermo(interval=max(n_steps // 5, 100))
    sim.run(n_steps)

    # ── Reduce trajectory → CellResult ─────────────────────────────
    if len(trajectory) < 2:
        # Degenerate cell (insertion failed): record NaNs.
        return CellResult(
            G_s_soft=G_s_soft, G_s_stiff=G_s_stiff,
            Ca_ratio=G_s_stiff / G_s_soft,
            theta_soft=np.nan, theta_stiff=np.nan, d_theta=np.nan,
            mean_dx_soft=np.nan, mean_dy_soft=np.nan,
            mean_dx_stiff=np.nan, mean_dy_stiff=np.nan,
            lane_order=np.nan,
            Z00=np.nan, Z01=np.nan, Z10=np.nan, Z11=np.nan,
            n_soft_placed=0, n_stiff_placed=0,
            walltime_s=time.time() - t_start,
        )

    t0 = trajectory[0]
    tN = trajectory[-1]
    pos0  = t0["positions"]
    posN  = tN["positions"]
    types = t0["types"]
    dx    = posN[:, 0] - pos0[:, 0]
    dy    = posN[:, 1] - pos0[:, 1]
    dx    = np.where(dx < 0, dx + NX, dx)      # un-wrap periodic-x

    soft_mask  = types == 0
    stiff_mask = types == 1

    theta_soft  = float(np.degrees(np.arctan2(
        dy[soft_mask].mean(), dx[soft_mask].mean())))
    theta_stiff = float(np.degrees(np.arctan2(
        dy[stiff_mask].mean(), dx[stiff_mask].mean())))

    snap_final = SimulationSnapshot.from_simulation(sim.core)
    lo         = lane_order(snap_final, axis="x")
    pts        = per_type_contact_stats(snap_final, contact_cutoff=0.5)
    Z          = pts.Z_matrix
    n_per_type = pts.n_particles

    # Pack a per-cell history.npz so the analysis script can re-plot
    # any cell on demand without re-running the simulation.
    traj_xy = np.stack([t["positions"] for t in trajectory])
    np.savez(os.path.join(output_dir, "history.npz"),
              traj_xy=traj_xy, types=types,
              dx_per_p=dx, dy_per_p=dy,
              G_s_soft=G_s_soft, G_s_stiff=G_s_stiff,
              theta_soft=theta_soft, theta_stiff=theta_stiff,
              lane_order=float(lo) if not np.isnan(lo) else 0.0,
              Z_matrix=Z, n_per_type=n_per_type)

    return CellResult(
        G_s_soft=G_s_soft, G_s_stiff=G_s_stiff,
        Ca_ratio=G_s_stiff / G_s_soft,
        theta_soft=theta_soft, theta_stiff=theta_stiff,
        d_theta=abs(theta_soft - theta_stiff),
        mean_dx_soft=float(dx[soft_mask].mean()),
        mean_dy_soft=float(dy[soft_mask].mean()),
        mean_dx_stiff=float(dx[stiff_mask].mean()),
        mean_dy_stiff=float(dy[stiff_mask].mean()),
        lane_order=float(lo) if not np.isnan(lo) else 0.0,
        Z00=float(Z[0, 0]) if Z.shape == (2, 2) else 0.0,
        Z01=float(Z[0, 1]) if Z.shape == (2, 2) else 0.0,
        Z10=float(Z[1, 0]) if Z.shape == (2, 2) else 0.0,
        Z11=float(Z[1, 1]) if Z.shape == (2, 2) else 0.0,
        n_soft_placed=int(n_per_type[0]) if n_per_type.size > 0 else 0,
        n_stiff_placed=int(n_per_type[1]) if n_per_type.size > 1 else 0,
        walltime_s=time.time() - t_start,
    )


# =================================================================
# Sweep driver
# =================================================================


def _build_pairs(values: list[float], full_grid: bool) -> list[tuple[float, float]]:
    """Upper-triangle pairs (soft < stiff) by default; full grid when --full.

    The diagonal (soft = stiff) is kept either way — it's the
    no-contrast control where Δθ should be ~0.
    """
    pairs: list[tuple[float, float]] = []
    for i, gs_soft in enumerate(values):
        for j, gs_stiff in enumerate(values):
            if full_grid:
                pairs.append((gs_soft, gs_stiff))
            else:
                if i <= j:        # upper triangle + diagonal
                    pairs.append((gs_soft, gs_stiff))
    return pairs


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--smoke", action="store_true",
                    help="3-value grid × 1000 steps (~3 min total)")
    p.add_argument("--full", action="store_true",
                    help="5×5 grid (25 cells) instead of upper-triangle (15)")
    p.add_argument("--steps", type=int, default=None,
                    help="override per-cell N_STEPS (defaults: smoke=1000, prod=20000)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv if argv is not None else sys.argv[1:])
    values = SWEEP_VALUES_SMOKE if args.smoke else SWEEP_VALUES_FULL
    n_steps = args.steps if args.steps is not None else (1000 if args.smoke else 20000)
    pairs = _build_pairs(values, full_grid=args.full)

    sweep_root = script_dir / "sweep_out"
    sweep_root.mkdir(exist_ok=True)

    print(f"\nSweep: {len(pairs)} cells, n_steps={n_steps}\n")
    print(f"  G_s values: {values}")
    print(f"  ETA: ~{len(pairs) * (1.5 if args.smoke else 12):.0f} min "
           f"(rough estimate)\n")

    results: list[CellResult] = []
    for k, (gs_soft, gs_stiff) in enumerate(pairs):
        cell_dir = sweep_root / f"G{gs_soft:.3f}_S{gs_stiff:.3f}"
        cell_dir.mkdir(exist_ok=True)
        print(f"[{k + 1:2d}/{len(pairs)}] G_s_soft={gs_soft:.3f} "
               f"G_s_stiff={gs_stiff:.3f}  → {cell_dir.name}", flush=True)
        try:
            r = run_one_cell(G_s_soft=gs_soft, G_s_stiff=gs_stiff,
                              output_dir=str(cell_dir),
                              n_steps=n_steps)
            results.append(r)
            print(f"     Δθ = {r.d_theta:6.2f}°,  "
                   f"θ_soft = {r.theta_soft:+6.2f}°,  "
                   f"θ_stiff = {r.theta_stiff:+6.2f}°,  "
                   f"Φ_lane = {r.lane_order:+5.2f},  "
                   f"wall {r.walltime_s:.0f} s")
        except Exception as e:
            print(f"     FAILED: {e}")
            continue

    # ── Aggregate into sweep_results.npz ───────────────────────────
    if not results:
        print("No successful cells. Bailing.")
        return 1

    field_names = list(CellResult.__dataclass_fields__.keys())
    cols = {name: np.asarray([getattr(r, name) for r in results])
             for name in field_names}
    np.savez(script_dir / "sweep_results.npz", **cols,
              sweep_values=np.asarray(values),
              n_cells=len(results))

    print(f"\nWrote {script_dir / 'sweep_results.npz'} ({len(results)} cells)")
    print("Run sweep_analyse.py next to produce the heatmap.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
