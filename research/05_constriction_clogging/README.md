# Project 5 — Clogging of deformable capsules at a microfluidic constriction

A dense suspension of soft (Skalak) capsules is driven through an abrupt
rectangular neck. Depending on the aperture ratio **D/d** and the capsule
deformability (capillary number **Ca**), the capsules either pass continuously
or build a stable arch that **clogs** the neck.

This is the everyday-life jamming problem — grain in a hopper, a crowd at a
doorway, coffee grounds bridging a filter, hair + grease in a drain — in its
cleanest microfluidic form. It is also a **validation target with a published
answer key** (see below), and it exercises the part of the engine that
projects 03 (DLD) and 04 (stenosis) did not: collective contact mechanics and
arching.

## Why this project

1. **Validation.** Reproduce a known clog/no-clog phase boundary.
2. **New everyday-life result.** Layer the jamming diagnostics
   (`pysoftflow.analysis.jamming`) — contact number Z, force-network
   percolation — onto the deformable-capsule clogging problem, which the
   existing (mostly rigid-particle) clogging literature does not report.
3. **GNN data source (Flavor B).** Unlike the dilute stenosis carriers, this
   is a *dense, contacting* suspension, so each frame produces a contact-rich
   graph. The sweep doubles as a labelled training set for the
   graph-neural-network clog predictor (`pysoftflow.ml.build_graph_dataset`).

## Control parameters

| Symbol | Meaning | Set by |
|---|---|---|
| `D/d`  | aperture / capsule-diameter ratio | neck width (`--aperture`) |
| `Ca`   | capillary number `= ρ ν r γ̇ / G_s`, `γ̇ = 4 u_max / W` | Skalak shear modulus `G_s` (`--ca`) |

## Validation answer key (literature)

| Source | Result to reproduce | Dimension |
|---|---|---|
| Bielinski, Aouane, Harting & Kaoui, arXiv:2110.13299 (2021) | rigid (Ca→0) clog at **D/d ≈ 3**; soft (**Ca > 0.005**) no clog even at D/d→1; `T_evac ∝ (d/D)^0.827, Ca^-0.108` | 3-D LBM-IBM |
| Hong, Kohne, Morrell, Wang & Weeks, PRE 96, 062605 (2017) | soft-particle clogging statistics | **2-D** (dimension-matched) |
| Marin et al., PRE 97 (2018); Souzy, Zuriguel & Marin, PRE 101 (2020) | rigid clog threshold D/d ≈ 2–4 | experiment |
| Marin & Souzy, Annu. Rev. Fluid Mech. 57, 89 (2025) | sieving / bridging / deposition regimes | review |

**Honest caveat:** Bielinski is 3-D; 2-D arches form more easily, so the
absolute threshold here will likely sit nearer D/d ≈ 2. The pass criterion is
reproducing the *trend* and the *deformability-driven shift*, anchored to the
2-D Hong (2017) experiment — not the exact 3-D number.

## Novelty hook (beyond validation)

The "obstacle above the outlet reduces clogging by up to ×100" result
(*Comm. Phys.* 4, 2021) is established only for dry grains and crowds. Repeating
it for **deformable capsules in viscous flow** is open — the everyday
"why a pillar before a doorway helps a crowd escape" question.

## Files

| File | Purpose |
|---|---|
| `run.py`     | single (D/d, Ca) configuration |
| `sweep.py`   | the (D/d, Ca) phase-diagram grid + per-cell GNN graph dump |
| `analyse.py` | clog detection, evacuation time, neck contact number Z; builds the GNN graph dataset. `--self-test` validates the logic without the engine |
| `references.md` | bibliography (verify DOIs via Crossref before any manuscript) |

## Build & run

The C++ core must be built first (it is not committed):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DSOFTFLOW_OPENMP=ON
cmake --build build -j
pip install -e ".[dev]"
```

Then:

```bash
# single run
python research/05_constriction_clogging/run.py --aperture 2.0 --ca 0.01 --smoke

# full phase diagram
python research/05_constriction_clogging/sweep.py            # ~hours
python research/05_constriction_clogging/sweep.py --smoke    # quick check

# analysis (works on any particle_data.csv)
python research/05_constriction_clogging/analyse.py --csv vtk_clog/particle_data.csv --neck-x 185

# verify the analysis logic with NO engine required:
python research/05_constriction_clogging/analyse.py --self-test
```

## Calibration notes (do this on the first smoke runs)

The flow parameters in `run.py` are first-pass estimates from the Bielinski
2021 scaling:

- `u_max` / `BODY_FX` are set for Stokes-regime microfluidics (Re ≈ 0.1).
  Confirm the empty-channel centre-line velocity after the first run and adjust
  `BODY_FX` if the measured `u_max` differs.
- `capillary_to_shear_modulus()` maps Ca → `G_s`. Check that the softest cells
  (large Ca) deform without numerical blow-up and the stiffest (small Ca)
  behave as near-rigid.
- The neck length `NECK_X` and chamber `PLUG_X` may need widening so the
  upstream plug holds all `N_CELLS` capsules at the chosen packing.

## Feeding the GNN (Flavor B)

`analyse.py --save-graphs` and `sweep.py` write per-cell `graphs.npz` (one
graph per frame, label = clog flag) via `pysoftflow.ml.build_graph_dataset`.
Nodes = capsules (position, velocity, wall distance, type, radius, contact
degree); edges = near-contacts (gap, direction, relative speed, approach
rate). Load with `pysoftflow.ml.load_graph_dataset`; bridge to PyTorch with
`ParticleGraph.to_pyg()` once `torch_geometric` is installed.
