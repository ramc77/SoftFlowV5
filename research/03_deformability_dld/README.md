# Project 3 — Deformability-based sorting in a DLD pillar array

## Research question

> **In a deterministic-lateral-displacement (DLD) pillar array, does
> the critical separation size depend on capsule stiffness contrast?
> Can soft particles be sorted by `G_s` alone with *no size
> difference*?**

The classical DLD result (Huang 2004, Holm 2011) is purely
geometric: the critical bin depends only on the pillar diameter,
gap, and lateral shift. For **rigid** spheres, that's been verified.
The open question is whether **deformability alone**, with identical
particle size, is enough to drive sorting in the same geometry.

## Why publishable now

| Reference | Relevance |
|---|---|
| Huang et al., *Science* **304**, 987 (2004) | Original DLD invention |
| Holm et al., *Lab Chip* **11**, 1326 (2011) | Critical-size formula for rigid spheres |
| D'Avino & Maffettone, *J Fluid Mech* **782**, 213 (2015) | Review of softness effects |
| Henon et al., *Biomicrofluidics* **11**, 064108 (2017) | Softness shifts the critical bin |
| Vahidkhah & Bagchi, *Soft Matter* **11**, 2097 (2015) | Capsule shape effects in 3D |

The gap: a **systematic 2D parameter sweep** of
`(G_s_soft, G_s_stiff)` at **identical particle size** mapping to
DLD displacement angle θ_DLD and lane-order parameter Φ_lane.

## What's measured

| Quantity | Source | Headline relevance |
|---|---|---|
| `theta_soft`, `theta_stiff` | `atan(<dy>/<dx>)` per type from trajectory log | DLD angle per species |
| `Δθ = |theta_soft − theta_stiff|` | Derived | **Sorting power** — the headline |
| `lane_order` | Phase-3 `lane_order(snap, axis="x")` | Cross-stream segregation |
| `Z_matrix` | Phase-3 `per_type_contact_stats` | Inter-type vs intra-type contacts |
| Trajectory `(x, y)` per particle per timestep | Step-callback logger | Figure 1 |

## Design

- **Identical particle size** (`R_PART = 3.0` for both soft and stiff).
  Mass + drag scale identically; only membrane stiffness differs.
- **10× stiffness contrast** (`G_s_soft = 0.03`, `G_s_stiff = 0.30`).
- **DLD geometry**: 8×4 pillar array on a parallelogram lattice,
  pillar radius 4, column spacing 35, row spacing 14, row-to-row
  y-shift = 3.5 (this is the "DLD step").
- Soft + stiff capsules seeded **mixed** upstream of the array;
  trajectories logged every 25 steps; final state compared.

## Methods (for the paper)

- **D2Q9 LBM** with **interpolated bounce-back** on the curved
  pillar surfaces (Bouzidi-Firdaouss-Lallemand 2001).
- **IBM** with 2 iterations of multi-direct forcing.
- **Skalak membrane** for both species. In the sweep, the two types
  differ **only in `G_s`** — `C_skalak`, `k_bend`, `k_area`, and
  `k_perimeter` are held identical so the diagonal cell
  (`G_s_soft = G_s_stiff`) is a true no-contrast control where
  Δθ ≡ 0 by symmetry. The standalone `run.py` uses a more realistic
  parameter contrast (RBC-like soft, tumour-like stiff).
- **Lubrication corrections** at close approach to the pillars.
- **Phase-3 `lane_order`** and **`per_type_contact_stats`** for
  quantitative segregation.
- **Randomized-order singleton seeding.** To avoid the order bias
  inherent in two-pass RSA (first species claims best positions),
  the sweep builds a list of `N_STIFF` "stiff" + `N_SOFT` "soft"
  labels, shuffles with the run's RNG seed, and places one capsule
  at a time. Each species sees the same expected placement density
  at every step; species labels become statistically interchangeable
  on the diagonal.

## Dimensionless numbers

| Number | Definition | Value | Regime |
|---|---|---|---|
| Re | `ū H / ν` | ~50 | Mild inertial regime |
| Ca_soft | `μ ū / G_s_soft` | ~0.2 | Significant deformation |
| Ca_stiff | `μ ū / G_s_stiff` | ~0.02 | Quasi-rigid |
| Ca ratio | 10 | | Stiffness contrast |
| `R / R_pillar` | 0.75 | | Particle-pillar size ratio |
| `λ_DLD` (row shift / column spacing) | 0.1 | | Sets the rigid critical bin |

For rigid spheres at this `λ_DLD`, Holm's formula gives a critical
size `D_c ≈ 4.0 / sqrt(λ_DLD) ≈ 12.6` (in pillar-gap units; our
particles are sub-critical). **Below the critical bin, rigid
spheres should *flow through* unsorted.** Any θ_DLD difference we
observe is therefore *softness-driven*.

## Limitations to declare

- **2D only.** Real DLD uses 3D micropillars; the third dimension
  modifies the local stress and may amplify or damp the deformability
  effect.
- **Single stiffness contrast.** A parameter sweep over `G_s_ratio`
  is the natural next step.
- **No diffusion or Brownian motion.** Sub-µm-particle DLD would
  need a thermal noise term, which the current engine doesn't
  include.
- **Pillar array is finite.** A real DLD has thousands of pillars;
  our 8×4 grid is for proof-of-concept.

## Running

```bash
python research/03_deformability_dld/run.py     # ~10 min
python research/03_deformability_dld/analyse.py # figures
```

Output:

- `vtk_dld/{fluid,particles}/*.pvd` — open in ParaView. Trajectories
  through the pillar array are striking.
- `history.npz` — full trajectory tensor, per-particle (dx, dy),
  per-type θ_DLD, lane order, Z_ij matrix.
- `fig1_trajectories.png` — overlaid trajectories soft (blue) vs
  stiff (red) showing path-divergence in the array.
- `fig2_dld_angle.png` — per-particle θ_DLD histogram, two coloured
  bars. If the two distributions are separated, sorting works.
- `fig3_summary_bars.png` — mean θ_DLD per type, lane order, and
  the Z_ij contact-matrix heatmap.

## Parameter sweep — the main result for the paper

A single `run.py` is the proof-of-concept; the actual paper rests on
a **`(G_s_soft, G_s_stiff)` sweep** showing that Δθ grows with the
stiffness ratio. The harness is `sweep.py`:

```bash
# Smoke (~10 min): 3-value grid × 1000 steps × 6 upper-triangle cells
python research/03_deformability_dld/sweep.py --smoke

# Production (~3 h): 5-value grid × 20 000 steps × 15 cells
python research/03_deformability_dld/sweep.py

# Full 5×5 grid (~5 h): use this for the heatmap appendix
python research/03_deformability_dld/sweep.py --full

# After the sweep, produce the heatmaps
python research/03_deformability_dld/sweep_analyse.py
```

**Sweep grid** (default `SWEEP_VALUES_FULL` in `sweep.py`):

```
G_s values: [0.03, 0.10, 0.20, 0.40, 0.80]
```

These all stay stable at `max_lattice_force = 0.04`. The 0.03 / 0.80
contrast gives `Ca_ratio = 27`, which is well above the threshold
where Henon et al. 2017 reported sorting in 3D.

**Output layout** after `sweep.py`:

```
research/03_deformability_dld/
├── sweep_out/
│   ├── G0.030_S0.030/             diagonal (no-contrast control)
│   │   ├── vtk_dld/...            full per-cell trajectory
│   │   ├── history.npz            cell-level time series
│   │   └── config/run_manifest.json   Phase-1 provenance
│   ├── G0.030_S0.100/
│   ├── ...
├── sweep_results.npz              aggregated grid (master file)
└── after `sweep_analyse.py`:
    ├── fig_sweep_dtheta.png                     Δθ_raw heatmap (G_s_soft × G_s_stiff)
    ├── fig_sweep_dtheta_corrected.png           Δθ_corr heatmap — THE HEADLINE
    ├── fig_sweep_dtheta_corrected_vs_Ca.png     Δθ_corr vs Ca-ratio collapse
    ├── fig_sweep_dtheta_vs_Ca.png               Δθ_raw vs Ca-ratio
    ├── fig_sweep_lane_order.png                 Φ_lane heatmap
    └── fig_sweep_summary.txt                    tabular summary for supplement
```

Every cell carries its own `run_manifest.json` (git SHA, compiler
flags, RNG seed, fully resolved params), so each grid point is
independently reviewer-traceable.

### Diagonal-subtracted Δθ — why and how

The diagonal (`G_s_soft = G_s_stiff`) is a physical no-contrast
control: both species are mechanically identical, so any nonzero
Δθ on the diagonal is a residual *artifact* (finite-N statistics,
small seeding-position fluctuations, finite simulation time). The
off-diagonal Δθ_raw is the *artifact + signal* sum.

`sweep_analyse.py` subtracts the per-row diagonal baseline:

```
Δθ_corr(soft, stiff) = Δθ_raw(soft, stiff)  −  Δθ_raw(soft, soft)
```

`Δθ_corr` is the pure deformability-contrast contribution. It is
what enters the headline figure (`fig_sweep_dtheta_corrected.png`)
and the Δθ vs Ca-ratio collapse plot. The raw heatmap is kept in
the supplement for transparency.

## Suggested paper structure

1. **Introduction** — DLD's rigid-sphere story; gap on soft particles
   at identical size.
2. **Methods** — LBM-IBM + IBB pillars + Skalak; geometry table.
3. **Results** —
   - Trajectories show systematic lateral drift differences (Fig 1).
   - θ_DLD distribution clearly bimodal (Fig 2).
   - Sweep `Ca_soft × Ca_stiff` → sorting fidelity heatmap.
4. **Discussion** — Mechanism: stiff particles experience larger
   lateral force per pillar interaction, soft particles deform and
   slip past. Comparison with Henon 2017 / Vahidkhah-Bagchi 2015.
5. **Conclusions** — DLD can sort at fixed particle size by
   stiffness alone, with sorting fidelity controlled by Ca ratio.

## Target journals

| Journal | IF (2024) | Fit |
|---|---|---|
| *Lab on a Chip* | ~7 | Primary — DLD is a central LoC method |
| *Microfluidics and Nanofluidics* | ~3.2 | Strong — direct fit |
| *Phys Rev Fluids* | ~3.0 | Strong — capsule mechanics in flow |
| *Soft Matter* | ~3.5 | Backup — coarse-grained capsule sorting |
