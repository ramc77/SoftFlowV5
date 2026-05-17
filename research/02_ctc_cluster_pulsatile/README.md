# Project 2 — CTC cluster stability under pulsatile shear with catch-slip bonds

> **Coarse-grained mechano-chemical proxy for circulating-cell
> aggregation under flow. NOT a validated cancer model.**

## Research question

> **Does pulsatile shear (heart-beat-like sinusoidal modulation of
> the driving body force) stabilise or fragment CTC clusters
> compared to steady shear at the same mean velocity, given the
> catch-slip nature of cell-cell adhesion bonds?**

The hypothesis is genuinely uncertain a priori:

- **Stabilising route.** Catch-slip bonds spend more time in their
  optimal force window when shear is pulsatile — they
  *strengthen* at intermediate force. Pulsatility could push the
  bond population into the catch regime cyclically, raising the
  mean lifetime relative to a constant-force regime.
- **Fragmenting route.** Peak systolic force may exceed the slip
  transition threshold and tear bonds faster than diastolic
  re-binding can recover, net-destabilising the cluster.

Neither has been directly tested in a 2D coarse-grained simulation
that resolves individual bonds AND tracks cluster size AND
characterises flow-rate-drop / spanning-cluster events.

## Why publishable now

| Reference | Relevance |
|---|---|
| Aceto et al., *Cell* **158**, 1110 (2014) | CTC clusters drive metastasis ~50× more than single cells |
| Au et al., *PNAS* **113**, 4947 (2016) | Clusters survive 5 % of normal circulation despite > critical shear |
| Thomas, Vogel, Sokurenko, *Annu Rev Biophys* **37**, 399 (2008) | Catch-slip bonds well-characterised under *constant* force only |
| Mol-mech of catch bonds, Pereverzev et al., *Biophys J* **89**, 1446 (2005) | Single-molecule kinetics |

The gap: the constant-force literature is rich, but **physiological
flow is pulsatile**. No published 2D coarse-grained model probes the
emergent question of whether pulsatility stabilises or fragments
multi-cell adhesive aggregates.

## What's measured

| Quantity | Source | Headline relevance |
|---|---|---|
| `largest_cluster` (steady, pulsatile) | Hoshen-Kopelman on the bond graph each `OUT_EVERY` steps | Aggregate persistence |
| `n_bonds` (steady, pulsatile) | `len(snap.bonds)` | Bond population dynamics |
| `Q_ratio` | EmbolizationDetector cross-section flow rate | Hemodynamic effect of the cluster |
| Peak vs mean cluster size | Time-series stats | Headline comparison |

## Design

The single script runs **two simulations back-to-back** with
identical seeds and identical *mean* body force:

- **Steady control**: `F_x(t) = BODY_FX_MEAN`
- **Pulsatile case**: `F_x(t) = BODY_FX_MEAN · (1 + amplitude · sin(2π t / period))`

with `amplitude = 0.6` (60 % modulation, comparable to physiological
systolic-diastolic ratios) and `period = 2000 steps` (a fraction of
one box-traversal time).

Outputs are written to **separate** `vtk_ctc_steady/` and
`vtk_ctc_pulsatile/` directories so they can be loaded side by side
in ParaView.

## Methods (for the paper)

- **2D D2Q9 LBM**, regularised collision (Latt & Chopard 2006).
- **Direct-forcing IBM** with 2 iterations of multi-direct forcing.
- **Skalak membrane** for all CTCs, `G_s = 0.10`, area conservation
  `C_skalak = 10`.
- **Bell + catch-slip adhesion** (Thomas et al. 2008):
  `k_off(F) = k_catch · exp(−F/F_catch) + k_slip · exp(F/F_slip)`
  with `k_catch = 0.05`, `F_catch = 0.015`, `k_slip = 0.001`,
  `F_slip = 0.012`.
- **Hoshen-Kopelman** on the bond graph for cluster identification.
- **Pulsatile body force** wired through `Simulation::setStepCallback`
  — the *only* C++ side-effect is the body-force update; the rest
  of the LBM-IBM stack is untouched.

## Dimensionless numbers

| Number | Definition | Value | Regime |
|---|---|---|---|
| Re | `ū H / ν` | ~30 | Inertial-viscous balance |
| Strouhal | `f · H / ū` | ~0.5 | Pulsation overlaps box-traversal time |
| Modulation depth | `(F_max − F_min) / 2 F_mean` | 0.6 | Strong physiological pulse |
| Ca | `μ ū / G_s` | ~0.05 | Modest deformation |

## Limitations to declare

- **Coarse-grained mechano-chemical proxy. Not a validated cancer
  model.** Repeated in every output and every figure caption.
- **2D only.** Real CTCs are 3D and the bond geometry is different.
- **Single CTC type.** Heterogeneity within real clusters (Aceto
  showed clusters are polyclonal) is not modelled.
- **No bond stochasticity beyond Poisson off-rate.** Real catch-slip
  bonds have richer two-state dynamics.
- **Pulsatile shape is sinusoidal.** Real arterial waveforms have a
  sharper systolic peak — easy to swap in once the sinusoidal case
  is established.

## Running

```bash
python research/02_ctc_cluster_pulsatile/run.py    # ~5-10 min
python research/02_ctc_cluster_pulsatile/analyse.py
```

Output:

- `vtk_ctc_steady/{fluid,particles}/*.pvd` — load steady in ParaView.
- `vtk_ctc_pulsatile/{fluid,particles}/*.pvd` — load pulsatile.
- `history.npz` — `largest_cluster`, `n_bonds`, `Q_ratio`,
  `body_fx`, `cluster_y_span` for both conditions.
- `fig1_cluster_size.png` — cluster size vs time (the headline).
- `fig2_bonds_vs_force.png` — bond count phase portrait.
- `fig3_summary_bars.png` — peak / mean comparison.

## Suggested paper structure

1. **Introduction** — CTC clusters drive metastasis; catch-slip
   bond literature; gap on pulsatile flow.
2. **Methods** — LBM-IBM + Skalak + Bell/catch-slip; V&V (Phase-1
   Poiseuille).
3. **Results** —
   - Cluster size vs time, steady vs pulsatile (Fig 1).
   - Bond population phase portrait (Fig 2): does pulsatility
     concentrate residence time in the catch-regime force window?
   - Peak / mean / fragmentation-event comparison (Fig 3).
   - Sweep: amplitude × period grid (production run).
4. **Discussion** — Mechanism: which catch/slip pathway dominates,
   how it depends on pulse parameters, biological implications.
5. **Conclusions** — Pulsatility shifts the steady-state bond
   distribution, with direction depending on whether peak force
   stays inside the catch window.

## Target journals

| Journal | IF (2024) | Fit |
|---|---|---|
| *Biophysical Journal* | ~3.6 | Strong — bond kinetics + cell mechanics |
| *Phys Rev Fluids* | ~3.0 | Strong — pulsatile flow + soft particle suspension |
| *PNAS* | ~10 | Possible if framed as "metastasis mechanics" |
| *Soft Matter* | ~3.5 | Backup — coarse-grained adhesion under flow |
