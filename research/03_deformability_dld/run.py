"""
Deformability-based sorting in a DLD-style pillar array
=======================================================
Research question
-----------------
In a deterministic-lateral-displacement (DLD) pillar array does the
critical separation size depend on capsule **stiffness contrast**?
Can soft particles be sorted by ``G_s`` alone with *no size
difference*?

Why this is publishable (2020-2025)
-----------------------------------
- Holm et al., *Lab Chip* **11**, 1326 (2011) derived the canonical
  DLD critical-size formula assuming rigid spheres.
- D'Avino & Maffettone, *J Fluid Mech* **782**, 213 (2015) reviewed
  soft-particle effects in microfluidic sorting.
- Henon et al., *Biomicrofluidics* **11**, 064108 (2017) and
  Vahidkhah & Bagchi, *Soft Matter* **11**, 2097 (2015) showed
  *softness* shifts the critical bin — but a systematic 2D study
  mapping (G_s_soft, G_s_stiff) → sorting fidelity is missing.

What's measured
---------------
- Lateral trajectory <y(x)> for each particle type — does the soft
  population displace less than the stiff one?
- Type-resolved exit-y distribution at the end of the pillar array.
- Lane order parameter Φ_lane (Phase 3) — quantifies the segregation
  *along* the array.
- Per-type contact matrix Z_ij at end of run.
- DLD displacement angle θ = atan(<Δy> / Δx) per type.

Output
------
- vtk_dld/{fluid,particles}/*.{vti,vtp}
- history.npz                     time series + per-particle final state
- config/run_manifest.json         provenance

References
----------
See ./references.md for the full citation list.
"""

import os
import sys

import numpy as np

# -- Path setup (find the build) --
script_dir  = os.path.dirname(os.path.abspath(__file__))
project_dir = os.path.dirname(os.path.dirname(script_dir))
sys.path.insert(0, os.path.join(project_dir, "python"))
for build_dir in ("build", "build_phase1", "build_phase2"):
    cand = os.path.join(project_dir, build_dir, "python")
    if os.path.isdir(cand):
        sys.path.insert(0, cand)

from pysoftflow import SoftFlowSimulation
from pysoftflow.analysis import SimulationSnapshot


# =================================================================
# Physical parameters (lattice units)
# =================================================================

NX, NY            = 400, 80
N_STEPS           = 800            # SMOKE; raise to 20 000 → 200 000 for production
TAU               = 0.7
BODY_FX           = 3e-5
MAX_LATTICE_FORCE = 0.04

# DLD pillar array: rows of circles, each row offset by `ROW_SHIFT`
# in y. The classical DLD geometry has pillars on a parallelogram
# lattice with period = N rows shifting cleanly back to the start.
N_PILLAR_COLS  = 8
N_PILLAR_ROWS  = 4
PILLAR_R       = 4.0
PILLAR_DX      = 35.0              # column spacing
PILLAR_DY      = 14.0              # row spacing
ROW_SHIFT      = 3.5               # row-to-row y-shift (the DLD step)
ARRAY_X_START  = 50.0
ARRAY_Y_CENTRE = NY / 2.0

# Bidisperse capsules — SAME size, different stiffness.
N_SOFT   = 20
N_STIFF  = 20
R_PART   = 3.0                     # IDENTICAL radius for both species

# Output cadence.
OUT_EVERY = 25
WARMUP_STEPS = 1000
WARMUP_RAMP  = 3000


# =================================================================
# Simulation setup
# =================================================================

sim = SoftFlowSimulation()

# -- Domain --
sim.domain(nx=NX, ny=NY)
sim.boundary(x="periodic", y="wall")
sim.fluid(tau=TAU, max_lattice_force=MAX_LATTICE_FORCE,
           collision="regularized")
sim.body_force(BODY_FX, 0.0)

# -- DLD pillar array: 8×4 grid of circles on a parallelogram lattice --
for col in range(N_PILLAR_COLS):
    for row in range(N_PILLAR_ROWS):
        cx = ARRAY_X_START + col * PILLAR_DX
        cy = (ARRAY_Y_CENTRE
              + (row - (N_PILLAR_ROWS - 1) / 2.0) * PILLAR_DY
              + col * ROW_SHIFT)
        # Keep pillars inside the channel.
        if cy < 6 or cy > NY - 6:
            continue
        sim.obstacle("circle", center=(cx, cy), radius=PILLAR_R)

sim.ibm(iterations=2)
sim.lubrication(enabled=True)

sim.output(format="vtk",
           directory=os.path.join(script_dir, "vtk_dld"),
           interval=OUT_EVERY)


# =================================================================
# Particle types: SAME size, different stiffness
# =================================================================

# Type 0: SOFT capsule.
sim.particle_type("soft",
                  model="skalak",
                  G_s=0.03,                  # 10× softer than the stiff
                  C_skalak=10.0,
                  k_bend=0.003, k_area=0.5, k_perimeter=0.05)

# Type 1: STIFF capsule (same size as soft, just stiffer membrane).
sim.particle_type("stiff",
                  model="skalak",
                  G_s=0.30,                  # 10× stiffer
                  C_skalak=10.0,
                  k_bend=0.020, k_area=1.0, k_perimeter=0.10)


# =================================================================
# Particle placement
# =================================================================

# Seed both species UPSTREAM of the pillar array, mixed in y.
sim.region("inlet",
            x=(10, ARRAY_X_START - 5),
            y=(8, NY - 8))

# Soft species first (random RSA fill), then stiff.
sim.generate("soft",  count=N_SOFT, region="inlet",
              radius=R_PART, num_nodes=18,
              method="random", seed=42,  min_gap=1.5)
sim.generate("stiff", count=N_STIFF, region="inlet",
              radius=R_PART, num_nodes=18,
              method="random", seed=43,  min_gap=1.5)


# =================================================================
# Trajectory logging via per-step callback
# =================================================================
#
# We sample (x, y) of every capsule's centroid every `OUT_EVERY`
# steps so analyse.py can plot trajectories and compute the DLD
# displacement angle per type.

trajectory: list[dict] = []          # one entry per snapshot


def _log_trajectory(core, step):
    if step % OUT_EVERY != 0 and step != 0:
        return
    snap = SimulationSnapshot.from_simulation(core)
    trajectory.append({
        "step": int(step),
        "time": float(step),
        "positions": snap.positions.copy(),
        "types":     snap.types.copy(),
        "radii":     snap.radii.copy(),
    })


sim.initialize()
sim.core.setStepCallback(_log_trajectory)

sim.warmup(steps=WARMUP_STEPS, ramp_steps=WARMUP_RAMP)
sim.thermo(interval=2000)
sim.run(N_STEPS)


# =================================================================
# Compute DLD displacement angle per type + save
# =================================================================

if len(trajectory) < 2:
    raise RuntimeError("not enough trajectory snapshots")

t0 = trajectory[0]
tN = trajectory[-1]

# Match capsule by index (the engine keeps capsules in registration
# order). Compute per-particle displacement and DLD angle.
positions_0 = t0["positions"]
positions_N = tN["positions"]
types_      = t0["types"]
# Periodic-x: we count *net* downstream travel including wraps. We
# don't unwrap explicitly because for a flow-driven run particles
# move only downstream; if a particle wraps, it just means it has
# done one full transit. Plot lateral displacement against the
# WRAPPED x for now.

dx_per_p = positions_N[:, 0] - positions_0[:, 0]
dy_per_p = positions_N[:, 1] - positions_0[:, 1]

# Account for periodic wrap in x: any particle that has wrapped will
# have negative dx_per_p (because final_x < initial_x); add NX to
# bring it into the positive range.
dx_per_p = np.where(dx_per_p < 0, dx_per_p + NX, dx_per_p)

soft_mask  = (types_ == 0)
stiff_mask = (types_ == 1)

theta_soft  = float(np.degrees(np.arctan2(dy_per_p[soft_mask].mean(),
                                              dx_per_p[soft_mask].mean())))
theta_stiff = float(np.degrees(np.arctan2(dy_per_p[stiff_mask].mean(),
                                              dx_per_p[stiff_mask].mean())))

# Phase-3 lane order at the end of run.
from pysoftflow.analysis.patterns import lane_order
snap_final = SimulationSnapshot.from_simulation(sim.core)
phi_lane   = lane_order(snap_final, axis="x")

# Phase-3 per-type contact matrix.
from pysoftflow.analysis.jamming import per_type_contact_stats
pts = per_type_contact_stats(snap_final, contact_cutoff=0.5)

# Pack trajectory into a 3-D ndarray (n_snapshots, n_particles, 2).
n_snap   = len(trajectory)
n_part   = positions_0.shape[0]
traj_xy  = np.stack([t["positions"] for t in trajectory])
traj_t   = np.array([t["time"]      for t in trajectory], dtype=np.float64)

np.savez(os.path.join(script_dir, "history.npz"),
          traj_xy=traj_xy,
          traj_t=traj_t,
          types=types_,
          radii=t0["radii"],
          dx_per_p=dx_per_p,
          dy_per_p=dy_per_p,
          theta_soft=theta_soft,
          theta_stiff=theta_stiff,
          lane_order=(float(phi_lane) if not np.isnan(phi_lane) else 0.0),
          Z_matrix=pts.Z_matrix,
          n_particles_per_type=pts.n_particles)

print("\n=== Headline numbers ===")
print(f"  soft  particles: <dx> = {dx_per_p[soft_mask].mean():.1f}, "
      f"<dy> = {dy_per_p[soft_mask].mean():+.2f}, "
      f"θ_DLD = {theta_soft:+.2f}°")
print(f"  stiff particles: <dx> = {dx_per_p[stiff_mask].mean():.1f}, "
      f"<dy> = {dy_per_p[stiff_mask].mean():+.2f}, "
      f"θ_DLD = {theta_stiff:+.2f}°")
print(f"  Δθ between species:    {abs(theta_soft - theta_stiff):.2f}°")
print(f"  lane-order parameter:   {phi_lane:+.3f}")
print(f"  per-type contact matrix Z_ij:")
print(f"    Z_soft-soft  = {pts.Z_matrix[0, 0]:.2f}")
print(f"    Z_soft-stiff = {pts.Z_matrix[0, 1]:.2f}")
print(f"    Z_stiff-soft = {pts.Z_matrix[1, 0]:.2f}")
print(f"    Z_stiff-stiff= {pts.Z_matrix[1, 1]:.2f}")
print(f"  → If Δθ > ~2° at the same particle size, DLD has sorted "
      f"by stiffness alone.")
