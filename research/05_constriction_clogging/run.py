"""
Project 5 — Clogging of deformable capsules at a microfluidic constriction
===========================================================================
A dense suspension of soft Skalak capsules is driven through an abrupt
rectangular constriction of aperture ``D``. Depending on the
aperture-to-particle ratio ``D/d`` and the capsule deformability
(capillary number ``Ca``), the capsules either pass continuously or
build a stable arch that *clogs* the neck.

This is the everyday-life "jamming at a doorway / grain in a hopper /
coffee-filter bridging" problem, and a validation target with a known
answer key:

  • Rigid limit (Ca -> 0):   permanent clog for  D/d  <~ 3
  • Soft  (Ca > ~0.005):     no clog even at  D/d -> 1
                              (Bielinski, Aouane, Harting & Kaoui,
                               arXiv:2110.13299, 2021 — 3-D LBM-IBM)
  • 2-D experimental anchor: Hong, Kohne, Morrell, Wang & Weeks,
                              Phys. Rev. E 96, 062605 (2017)
  • Constricted-suspension experiment: Marin et al., PRE 97 (2018);
                              Souzy, Zuriguel & Marin, PRE 101 (2020)

Control parameters
------------------
  D/d   — aperture ratio   (set by the neck width; --aperture)
  Ca    — capillary number (set by the Skalak shear modulus; --ca)
          Ca = rho * nu * r * gamma_dot / G_s,  gamma_dot = 4 u_max / W

Output (in ./vtk_clog/):
  • particle_data.csv          per-step capsule centroids + velocities
                               (consumed by analyse.py and the GNN
                                graph extractor pysoftflow.ml.build_graph_dataset)
  • vtk_clog/{fluid,particles} ParaView
  • config/run_manifest.json   provenance

NOTE: parameters below are first-pass estimates from the Bielinski 2021
scaling; calibrate u_max / G_s on the first smoke runs (see README).

Usage:
    python research/05_constriction_clogging/run.py --aperture 2.0 --ca 0.01
    python research/05_constriction_clogging/run.py --smoke
"""

from __future__ import annotations

import argparse
import pathlib
import sys

import numpy as np

# -- Path setup (mirror project 04) --
script_dir = pathlib.Path(__file__).resolve().parent
project_dir = script_dir.parents[1]
sys.path.insert(0, str(project_dir / "python"))
for d in ("build", "build_phase1", "build_phase2"):
    cand = project_dir / d / "python"
    if cand.is_dir():
        sys.path.insert(0, str(cand))

from pysoftflow import SoftFlowSimulation  # noqa: E402

# =================================================================
# Physical parameters (lattice units)
# =================================================================

NX, NY = 220, 48          # channel: length x width (W = NY)
# A narrower channel + higher viscosity shortens the diffusive flow-
# development time tau_dev ~ W^2/nu (here ~48^2/0.133 ~ 1.7e4 steps), so a
# feasible smoke run reaches developed flow.
TAU = 0.9                 # nu = (tau - 0.5)/3 = 0.1333
NU = (TAU - 0.5) / 3.0

# Capsule suspension
D_CAPSULE = 8.0           # rest diameter d  (radius r = 4)
R_CAPSULE = D_CAPSULE / 2.0
N_CELLS = 40              # dense plug (generator places as many as fit)
NUM_NODES = 16            # membrane resolution per capsule
PLUG_X = (40, 145)        # plug seated AT the neck mouth so the leading
PLUG_Y = (8, NY - 8)      # capsules reach it almost immediately (centres >= 2r
                          # from the walls to avoid wall overlap)

# Constriction: two rectangular obstacles leaving a central gap D.
NECK_X = (150, 175)       # streamwise extent of the neck (centre = 162.5)
APERTURE_DEFAULT = 2.0    # D/d
CA_DEFAULT = 0.01         # capillary number

# Flow driving: PERIODIC channel + body force (the stable, well-tested path;
# this configuration ran 32k steps without divergence at D/d=2.0, producing a
# real clog with neck contacts). The velocity-driven inlet/outlet variant is
# more controlled across apertures but proved numerically fragile at the
# extreme necks, so for data generation we use the robust body-force drive.
# Gentle drive keeps Re ~ 0.24, Ma ~ 0.014 (stable envelope).
TARGET_U_MAX = 0.008
U_MAX = TARGET_U_MAX
_FX_CAL, _UMAX_CAL = 4.63e-6, 1.65e-3         # measured calibration point
BODY_FX = _FX_CAL * (U_MAX / _UMAX_CAL)

# Contact / near-field physics (dense suspension needs these).
MORSE_EPS = 1.0e-3
MORSE_SIGMA = 1.0
MORSE_RCUT = 2.0
# NOTE: the engine also supports a DEM-style dissipative/frictional contact
# (style="dem" with damping_normal / friction_coeff; see docs/theory/contact.md
# and tests/python/test_contact_friction_damping.py). It is left OFF here while
# the inlet/outlet flow is being stabilised, to keep this run on the simplest
# (pure-repulsion) contact path.

OUT_EVERY = 200


def capillary_to_shear_modulus(ca: float) -> float:
    """Map a target capillary number to the Skalak shear modulus G_s.

    Ca = rho * nu * r * gamma_dot / G_s, with rho = 1 and the wall shear
    rate gamma_dot = 4 u_max / W (Bielinski 2021 convention). Smaller Ca =>
    stiffer capsule (larger G_s).
    """
    gamma_dot = 4.0 * U_MAX / NY
    return NU * R_CAPSULE * gamma_dot / ca


def parse_args(argv):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--smoke", action="store_true", help="Short run (5000 steps)")
    p.add_argument("--aperture", type=float, default=APERTURE_DEFAULT,
                   help="Aperture ratio D/d (default 2.0)")
    p.add_argument("--ca", type=float, default=CA_DEFAULT,
                   help="Capillary number (default 0.01)")
    p.add_argument("--steps", type=int, default=None, help="Override total steps")
    p.add_argument("--out", type=str, default=str(script_dir / "vtk_clog"),
                   help="Output directory")
    return p.parse_args(argv)


def build_sim(*, aperture: float, ca: float, output_dir: str, rng_seed: int = 0x05C106,
              placement: str = "fill", jitter: float = 0.0):
    """Construct a constriction-clogging simulation for given (D/d, Ca).

    ``placement`` selects the seeding method: ``"fill"`` (dense hex packing,
    default) or ``"random"`` (seed-dependent RSA). For multiple distinct seeds
    of a dense packing, use ``placement="fill"`` with ``jitter>0`` — this keeps
    the dense contacts (unlike RSA, which spreads capsules out) while letting
    different ``rng_seed`` values produce different configurations.
    """
    neck_width = aperture * D_CAPSULE
    gap_lo = (NY - neck_width) / 2.0
    gap_hi = (NY + neck_width) / 2.0
    g_s = capillary_to_shear_modulus(ca)

    sim = SoftFlowSimulation()
    sim.domain(nx=NX, ny=NY)
    # Periodic channel driven by a body force: capsules recirculate naturally
    # (periodic fluid + capsules), so the finite suspension sustains its
    # concentration at the neck. This is the stable data-generation path.
    sim.boundary(x="periodic", y="wall")
    sim.fluid(tau=TAU, collision="regularized")
    sim.body_force(BODY_FX, 0.0)
    sim.ibm(iterations=2)
    sim.lubrication(enabled=True, h_threshold=1.5)
    sim.interaction("particle-particle", style="morse",
                    epsilon=MORSE_EPS, sigma=MORSE_SIGMA, r_cut=MORSE_RCUT)
    sim.interaction("particle-obstacle", style="morse",
                    epsilon=MORSE_EPS, sigma=MORSE_SIGMA, r_cut=MORSE_RCUT)

    # Constriction: bottom + top rectangular obstacles leaving a gap D.
    sim.obstacle("rect", p1=(NECK_X[0], 0.0), p2=(NECK_X[1], gap_lo))
    sim.obstacle("rect", p1=(NECK_X[0], gap_hi), p2=(NECK_X[1], float(NY)))

    sim.output(format="vtk", directory=output_dir, interval=OUT_EVERY)

    sim.particle_type("cell", model="skalak",
                      G_s=g_s, C_skalak=10.0,
                      k_bend=0.02, k_area=1.0, k_perimeter=0.1)
    sim.region("inlet_plug", x=PLUG_X, y=PLUG_Y)
    if placement == "fill":
        sim.generate("cell", count=N_CELLS, region="inlet_plug",
                     radius=R_CAPSULE, num_nodes=NUM_NODES,
                     method="fill", packing="hex", seed=rng_seed, min_gap=1.0,
                     jitter=jitter)
    else:
        sim.generate("cell", count=N_CELLS, region="inlet_plug",
                     radius=R_CAPSULE, num_nodes=NUM_NODES,
                     method=placement, seed=rng_seed, min_gap=1.0)

    sim.initialize()
    return sim, dict(aperture=aperture, ca=ca, neck_width=neck_width, g_s=g_s)


def main(argv=None) -> int:
    args = parse_args(argv if argv is not None else sys.argv[1:])
    n_steps = args.steps if args.steps is not None else (25_000 if args.smoke else 80_000)

    sim, meta = build_sim(aperture=args.aperture, ca=args.ca, output_dir=args.out)
    print("\nProject 5 — constriction clogging")
    print(f"  D/d        = {meta['aperture']:.2f}  (neck width = {meta['neck_width']:.1f} LU)")
    print(f"  Ca         = {meta['ca']:.4g}  (G_s = {meta['g_s']:.4g})")
    print(f"  u_max~{U_MAX:.4g} LU,  body_fx={BODY_FX:.3g},  "
          f"Re ~ {U_MAX * R_CAPSULE / NU:.2f}  (periodic + body force)")
    print(f"  N capsules = {N_CELLS},  steps = {n_steps},  out = {args.out}\n")

    sim.thermo(interval=max(n_steps // 25, 100))
    sim.warmup(steps=500, ramp_steps=1500)
    sim.run(n_steps)
    print("\nDone. Analyse with: python research/05_constriction_clogging/analyse.py "
          f"--csv {args.out}/particle_data.csv --neck-x {0.5 * (NECK_X[0] + NECK_X[1]):.0f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
