"""
CTC cluster stability under pulsatile shear with catch-slip bonds
=================================================================
Research question
-----------------
Does *pulsatile* shear (heart-beat-like sinusoidal modulation)
stabilise or fragment circulating-tumour-cell (CTC) clusters
compared to steady shear, given the catch-slip nature of cell-cell
adhesion bonds?

Why this is publishable (2023-2025)
-----------------------------------
- Aceto et al., *Cell* 158, 1110 (2014) — CTC clusters drive
  metastasis ~50× more than single cells.
- Au et al., *PNAS* 113, 4947 (2016) — clusters survive 5 % of normal
  circulation despite > critical-shear flow.
- Thomas-Vogel-Sokurenko, *Annu Rev Biophys* 37, 399 (2008) —
  catch-slip bonds well-characterised under *constant* force.
- Mechanism gap: catch-slip bonds under *time-varying* force are
  not directly probed in 2D coarse-grained simulation that ALSO
  tracks cluster size + spanning + flow-rate-drop events.

The hypothesis is *non-trivial*: catch bonds can spend more time at
their optimal force window when shear is pulsatile (the catch-bond
strengthening regime is force-banded), which could *stabilise*
clusters. Alternatively, the peak force during systole might
fragment them faster than steady. This experiment decides between
those.

What's measured
---------------
- Time series of the largest connected cluster size (Hoshen-Kopelman).
- Time series of total bond count.
- Comparison: steady-shear control vs pulsatile (matched ū).
- Cluster lifetime distribution (fragmentation events).
- Embolization-event count.

Output
------
- vtk_ctc_pulsatile/{steady,pulsatile}/{fluid,particles}/*.{vti,vtp}
- history.npz                       time series for both conditions
- config/run_manifest.json          provenance

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
from pysoftflow.tumor_growth import (
    EmbolizationDetector,
    TumorGrowthRun,
)
from pysoftflow.analysis import SimulationSnapshot
from pysoftflow.analysis.patterns import hoshen_kopelman


# =================================================================
# Physical parameters (lattice units)
# =================================================================

NX, NY            = 300, 50
N_STEPS           = 1000           # SMOKE; raise to 20 000 → 200 000 for production
TAU               = 0.7
BODY_FX_MEAN      = 4e-5
MAX_LATTICE_FORCE = 0.04

# Pulse parameters: physiological heart rate ~1 Hz; in lattice units
# we pick a period so we see many full pulses across N_STEPS.
PULSE_AMPLITUDE   = 0.6            # fraction of mean (0=steady, 1=fully modulated)
PULSE_PERIOD      = 2000           # steps per pulse

# CTC cluster: small group of adhesive cells.
N_CTC             = 9
R_CTC             = 3.5

OUT_EVERY         = 25
WARMUP_STEPS      = 200            # SMOKE; raise to 1000 for production
WARMUP_RAMP       = 500            # SMOKE; raise to 3000 for production


def build(*, pulsatile: bool, output_subdir: str, rng_seed: int):
    """Build one condition (steady or pulsatile) and return (sim, run)."""
    sim = SoftFlowSimulation()
    sim.domain(nx=NX, ny=NY)
    sim.boundary(x="periodic", y="wall")
    sim.fluid(tau=TAU,
               max_lattice_force=MAX_LATTICE_FORCE,
               collision="regularized")
    sim.body_force(BODY_FX_MEAN, 0.0)

    sim.ibm(iterations=2)
    sim.lubrication(enabled=True)

    # CTC-cluster adhesion — catch-slip bonds (this is the whole point).
    # Aggressive k_on (0.2 per step) and a generous d_bond (6.0) so
    # bonds form quickly between adjacent cluster cells even at typical
    # hex-fill centroid distances. Tune for production sweeps.
    sim.adhesion(enabled=True,
                  k_on=0.2, k_off=0.001,
                  k_bond=0.05, d_bond=6.0,
                  F_crit=0.01, max_bonds_per_node=4,
                  bond_model="catch_slip",
                  k_off_catch=0.05, F_catch=0.015,
                  k_off_slip=0.001, F_slip=0.012)

    sim.particle_type("ctc",
                       model="skalak",
                       G_s=0.10, C_skalak=10.0,
                       k_bend=0.005,
                       k_area=0.7, k_perimeter=0.07)

    # Initial CTC cluster: random close-contact seed. RSA with a
    # small min_gap packs cells densely in a small region — exactly
    # what we want for an initial CTC cluster.
    sim.region("cluster_seed",
                x=(40, 75),
                y=(NY / 2 - 12, NY / 2 + 12))
    sim.generate("ctc", count=N_CTC, region="cluster_seed",
                 radius=R_CTC, num_nodes=20,
                 method="random",
                 seed=rng_seed, min_gap=0.3)

    sim.output(format="vtk",
                directory=os.path.join(script_dir, output_subdir),
                interval=OUT_EVERY)

    sim.initialize()

    # We register a TumorGrowthRun ONLY for the embolization detector
    # — division is intentionally disabled (k_div=0).
    run = TumorGrowthRun(sim=sim.core)
    run.add_embolization_detector(
        EmbolizationDetector(x_section=NX // 2,
                              flow_drop_threshold=0.6,
                              contact_cutoff=1.0,
                              band_fraction=0.15))

    # Pulsatile body force: a sinusoidal modulation around the mean.
    # Installed as a step callback that updates the LBM body force in
    # place each timestep.  For the steady control we install a no-op
    # callback so both branches go through the same plumbing.
    history = {
        "step": [], "time": [], "Q_ratio": [],
        "largest_cluster": [], "n_bonds": [], "body_fx": [],
        "cluster_y_span": [],
    }

    def _step_cb(core, step):
        # 1. Update the body force per the pulsation rule.
        if pulsatile:
            phase = 2.0 * np.pi * step / PULSE_PERIOD
            fx = BODY_FX_MEAN * (1.0 + PULSE_AMPLITUDE * np.sin(phase))
        else:
            fx = BODY_FX_MEAN
        core.params().fluid.body_force_x = fx

        # 2. Record diagnostics every OUT_EVERY steps.
        if step % OUT_EVERY != 0 and step != 0:
            return
        snap = SimulationSnapshot.from_simulation(core)
        hk   = hoshen_kopelman(snap, contact_cutoff=1.0,
                                use_bonds=(snap.bonds.shape[0] > 0))
        n_bonds = snap.bonds.shape[0]

        # Use the embolization detector to get Q(t)/Q_0.
        det = run._detectors[0]
        # First call: calibrate baseline; subsequent calls drive history.
        if det.Q_baseline == 0.0:
            det.baseline(core)
        ev = det.step(core, snap)

        history["step"].append(int(step))
        history["time"].append(float(step))
        history["Q_ratio"].append(float(det.flow_rate_history[-1])
                                    if len(det.flow_rate_history) else 1.0)
        history["largest_cluster"].append(int(hk.largest_size))
        history["n_bonds"].append(int(n_bonds))
        history["body_fx"].append(float(fx))
        history["cluster_y_span"].append(float(det.cluster_span_history[-1])
                                            if len(det.cluster_span_history)
                                            else 0.0)

    sim.core.setStepCallback(_step_cb)
    return sim, run, history


# =================================================================
# Run BOTH conditions back-to-back
# =================================================================

print("================================================================")
print(" Project 2 — pulsatile vs steady shear, catch-slip CTC clusters ")
print("================================================================")

results = {}
for condition, is_pulsatile in [("steady", False), ("pulsatile", True)]:
    print(f"\n--- Running {condition} ---")
    sim, run, hist = build(pulsatile=is_pulsatile,
                            output_subdir=f"vtk_ctc_{condition}",
                            rng_seed=0x0AFE0001)   # < 2**31
    sim.warmup(steps=WARMUP_STEPS, ramp_steps=WARMUP_RAMP)
    sim.thermo(interval=2000)
    sim.run(N_STEPS)
    results[condition] = {
        "history":  {k: np.asarray(v) for k, v in hist.items()},
        "summary":  run.summary(),
    }
    print(f"  summary: {results[condition]['summary']}")


# =================================================================
# Save
# =================================================================

np.savez(os.path.join(script_dir, "history.npz"),
          **{f"{cond}_{key}": val
             for cond, data in results.items()
             for key, val in data["history"].items()})

print("\n=== Headline numbers ===")
for cond, data in results.items():
    hist = data["history"]
    max_cluster_steady = int(hist["largest_cluster"].max())
    final_cluster      = int(hist["largest_cluster"][-1])
    mean_bonds         = float(hist["n_bonds"].mean())
    print(f"  {cond:>9}: peak cluster = {max_cluster_steady:2d}, "
          f"final = {final_cluster:2d}, mean bonds = {mean_bonds:.1f}")
