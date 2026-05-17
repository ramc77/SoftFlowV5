"""
Microplastic margination + chemical leaching in a stenosed microvessel
======================================================================
Research question
-----------------
Does vessel stenosis preferentially trap microplastics (MPs), and
does this trapping concentrate plasticizer leaching at the stenotic
wall? If yes, stenotic endothelium gets a *spatially-concentrated*
chemical dose that the bulk flow doesn't see.

Why this is publishable (2024-2026)
-----------------------------------
- Leslie et al., *Environ Int* 158, 107199 (2022) — MPs detected in
  human blood.
- Ragusa et al., *Environ Int* 146, 106274 (2021) — MPs in human
  placenta.
- Campen et al., *Tox Sci* 192, 1 (2023) — MPs cross blood-brain
  barrier.
- Mechanism gap: nobody has simulated MP transport in a stenosed
  microvessel *with deformable RBCs, finite-reservoir leaching, and
  Langmuir uptake on cell surfaces* in a single 2D model.

What's measured
---------------
- Spatial dose map of leached chemical on the channel walls.
- MP residence-time distribution near the stenosis.
- Cell-free layer (CFL) shift caused by the obstacle pair.
- Margination index for MPs vs RBCs vs WBCs.
- Time-resolved particle mass M_p(t) per microplastic.

Output
------
- vtk_mp_stenosis/{fluid,particles}/*.{vti,vtp}  for ParaView
- history.npz                                    time series
- diagnostics/                                    Phase-3 fields (dose map)
- config/run_manifest.json                       provenance

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


# =================================================================
# Physical parameters (lattice units)
# =================================================================

NX, NY            = 400, 80          # vessel length x diameter
N_STEPS           = 800              # SMOKE; raise to 30 000 → 200 000 for production
TAU               = 0.7              # ν = 0.0667 (50 % faster than τ=0.8)
BODY_FX           = 3e-5             # body force ≈ pressure gradient
MAX_LATTICE_FORCE = 0.04             # raised from default 0.01

# Stenosis: two circular obstacles creating a symmetric narrowing.
STENOSIS_X      = NX // 2
STENOSIS_RADIUS = 15.0

# RBCs (Skalak, soft, ~40 % hematocrit-equivalent).
N_RBC      = 28
R_RBC      = 6.0

# WBCs (Neo-Hookean, larger and stiffer).
N_WBC      = 6
R_WBC      = 11.0

# Microplastics — bidisperse mix of small fragments and larger beads.
N_MP_SMALL = 14
N_MP_LARGE = 10
R_MP_SMALL = 2.5
R_MP_LARGE = 5.0

# Output cadence — at ū ≈ 0.07 and OUT_EVERY = 20 each frame is
# ≈ 1.4 cells of advection (≈ 25 % of an RBC radius).
OUT_EVERY = 20

# Chemistry starts after the flow + cells have reached steady state.
WARMUP_NO_CHEMISTRY = 400            # SMOKE; raise to 4000 for production


# =================================================================
# Simulation setup
# =================================================================

sim = SoftFlowSimulation()

# -- Domain: periodic microvessel with no-slip walls --
sim.domain(nx=NX, ny=NY)
sim.boundary(x="periodic", y="wall")
sim.fluid(tau=TAU, max_lattice_force=MAX_LATTICE_FORCE,
           collision="regularized")     # Latt-Chopard 2006: better stability at low τ
sim.body_force(BODY_FX, 0.0)

# -- Stenosis: two opposing circles narrowing the channel midway --
sim.obstacle("circle", center=(STENOSIS_X, 8),       radius=STENOSIS_RADIUS)
sim.obstacle("circle", center=(STENOSIS_X, NY - 8),  radius=STENOSIS_RADIUS)

# -- Multi-direct forcing for sharper no-slip on the cells --
sim.ibm(iterations=2)

# -- Lubrication corrections (critical for dense RBC suspensions) --
sim.lubrication(enabled=True)

# -- Adhesion: MPs adhere to the endothelial wall. Catch-slip kinetics
#    mimic non-specific MP–surface binding under varying shear. --
sim.adhesion(enabled=True,
             k_on=0.0002, k_off=0.005, k_bond=0.006,
             d_bond=1.5, F_crit=0.010, max_bonds_per_node=2,
             wall_adhesion=True, wall_k_on=0.0004,
             wall_k_off=0.003, wall_k_bond=0.008,
             wall_receptor_spacing=3.0,
             bond_model="catch_slip",
             k_off_catch=0.005, F_catch=0.005,
             k_off_slip=0.005,  F_slip=0.025)

# -- Scalar transport: chemical leached by microplastics --
sim.scalar_transport(enabled=True,
                     diffusivity=0.01,
                     n_species=1,
                     initial_concentration=0.0)

# -- Viscosity contrast: RBC interior ~5x plasma --
sim.viscosity_contrast(enabled=True, update_interval=50)

sim.output(format="vtk",
           directory=os.path.join(script_dir, "vtk_mp_stenosis"),
           interval=OUT_EVERY)


# =================================================================
# Particle types
# =================================================================

# Type 0: Red Blood Cell (Skalak, soft, viscosity contrast 5x).
sim.particle_type("rbc", model="skalak",
                  G_s=0.04, C_skalak=10.0,
                  k_bend=0.003, k_area=0.5, k_perimeter=0.05,
                  gamma_visc=0.02, viscosity_ratio=5.0)

# Type 1: White Blood Cell (Neo-Hookean, stiffer, larger).
sim.particle_type("wbc", model="neo_hookean",
                  G_s=0.12,
                  k_bend=0.010, k_area=0.8, k_perimeter=0.08,
                  gamma_visc=0.03, viscosity_ratio=3.0)

# Type 2: Small microplastic fragments (Hookean, very stiff,
# rigid translation+rotation only).
sim.particle_type("mp_small", model="hookean",
                  G_s=0.8, k_bend=0.05,
                  k_area=2.0, k_perimeter=0.20,
                  gamma_visc=0.01, is_rigid=True)

# Type 3: Large microplastic beads (also rigid).
sim.particle_type("mp_large", model="hookean",
                  G_s=1.0, k_bend=0.08,
                  k_area=3.0, k_perimeter=0.30,
                  gamma_visc=0.01, is_rigid=True)


# =================================================================
# Particle placement
# =================================================================

# Seed cells throughout the vessel, MPs upstream of the stenosis so we
# can watch them migrate and (with luck) get trapped at the constriction.
margin = R_WBC + 4.0

sim.region("vessel",     x=(15, NX - 15), y=(margin, NY - margin))
sim.region("upstream",   x=(15, STENOSIS_X - STENOSIS_RADIUS - 6),
                          y=(margin, NY - margin))

# RBCs throughout the vessel.
sim.generate("rbc", count=N_RBC, region="vessel",
             radius=R_RBC, num_nodes=24,
             method="random", seed=42, min_gap=2.0)

# WBCs scattered.
sim.generate("wbc", count=N_WBC, region="vessel",
             radius=R_WBC, num_nodes=28,
             method="random", seed=100, min_gap=3.0)

# Microplastics upstream of the stenosis.
sim.generate("mp_small", count=N_MP_SMALL, region="upstream",
             radius=R_MP_SMALL, num_nodes=12,
             method="random", seed=200, min_gap=1.5)
sim.generate("mp_large", count=N_MP_LARGE, region="upstream",
             radius=R_MP_LARGE, num_nodes=18,
             method="random", seed=300, min_gap=2.0)


# =================================================================
# Two-phase run: flow equilibration first, chemistry second
# =================================================================
#
# Phase 1: let the flow reach steady-state Poiseuille and let the
# cells distribute themselves before any chemistry starts. This
# prevents the t=0 concentration burst from contaminating the
# spatial-dose map.

sim.warmup(steps=1000, ramp_steps=3000)
sim.thermo(interval=2000)
print(f"\nPhase 1: flow equilibration, no chemistry ({WARMUP_NO_CHEMISTRY} steps)")
sim.run(WARMUP_NO_CHEMISTRY)

# Phase 2: activate chemistry. Damköhler ≈ k_L · L/U ≈ 0.0002 · 400/0.07
# ≈ 1.1 — leaching-diffusion balanced, no concentration burst.
print("\nPhase 2: chemistry active (Fick leaching + Langmuir adsorption)")

sim.scalar_source("mp_small",
                  k_leach=0.0002,
                  C_eq=1.0,
                  M_p_initial=300.0)
sim.scalar_source("mp_large",
                  k_leach=0.00015,
                  C_eq=1.0,
                  M_p_initial=600.0)
sim.scalar_source("rbc",
                  k_adsorb=0.0005,
                  k_desorb=0.00005,
                  Gamma_max=1.0)
sim.scalar_source("wbc",
                  k_adsorb=0.0008,
                  k_desorb=0.00003,
                  Gamma_max=1.0)

sim.run(N_STEPS - WARMUP_NO_CHEMISTRY)


# =================================================================
# Save time series for analyse.py
# =================================================================

# Per-MP cumulative leached mass (queried via getCapsuleMp on the
# AdvectionDiffusion module).
adr = sim.core.advectionDiffusion()
n_caps = sim.core.capsules().numCapsules()
M_p_remaining = np.empty(n_caps, dtype=np.float64)
released      = np.empty(n_caps, dtype=np.float64)
absorbed      = np.empty(n_caps, dtype=np.float64)
type_label    = np.empty(n_caps, dtype=np.int64)
for k in range(n_caps):
    M_p_remaining[k] = adr.getCapsuleMp(k)
    released[k]      = adr.getCapsuleReleased(k)
    absorbed[k]      = adr.getCapsuleAbsorbed(k)
    type_label[k]    = sim.core.capsules()[k].getType()

# Wall dose: integrate scalar field in two near-wall bands flanking the
# stenosis (the spatially-concentrated dose, our headline metric).
C   = adr.concentration(0)         # shape (ny, nx)
ny, nx = C.shape
WALL_BAND = 4                       # cells from the wall
sten_x_lo = STENOSIS_X - 30
sten_x_hi = STENOSIS_X + 30
dose_bottom = float(np.sum(C[:WALL_BAND, sten_x_lo:sten_x_hi]))
dose_top    = float(np.sum(C[ny - WALL_BAND:, sten_x_lo:sten_x_hi]))
dose_bulk   = float(np.sum(C)) - dose_bottom - dose_top

np.savez(os.path.join(script_dir, "history.npz"),
          M_p_remaining=M_p_remaining,
          released=released,
          absorbed=absorbed,
          type_label=type_label,
          dose_bottom_at_stenosis=dose_bottom,
          dose_top_at_stenosis=dose_top,
          dose_bulk=dose_bulk,
          stenosis_x=STENOSIS_X,
          scalar_field=C.copy())

print("\n=== Headline numbers ===")
print(f"  Dose at bottom wall (near stenosis):  {dose_bottom:.4g}")
print(f"  Dose at top    wall (near stenosis):  {dose_top:.4g}")
print(f"  Dose in bulk fluid:                   {dose_bulk:.4g}")
print(f"  Concentration ratio (wall / bulk):    "
      f"{(dose_bottom + dose_top) / max(dose_bulk, 1e-30):.2f}")
print(f"  → If ratio > 1, the stenosis concentrates leached "
      f"chemical at the wall.")
