# Project 1 — Microplastic margination + leaching in stenosed microvasculature

## Research question

> **Does vessel stenosis preferentially trap microplastics, and does
> this trapping spatially concentrate plasticizer leaching at the
> stenotic wall?**

A positive answer means that *atherosclerotic* / *constricted*
vessels — already disease-relevant — get a doubly-localised dose of
microplastic-borne chemicals (BPA, phthalates, PFAS), with two
clinically relevant consequences:

1. The endothelial cells near a stenosis see a chemical dose
   substantially higher than the bulk-flow concentration.
2. The MP residence time at the constriction is significantly
   prolonged vs unobstructed vessels, giving the leached chemicals
   more time to cross the endothelium.

## Why it's publishable now

Microplastics have been detected in human blood (Leslie et al.,
*Environ Int* 2022), in placenta (Ragusa et al., *Environ Int* 2021),
and in brain tissue (Campen et al., *Tox Sci* 2023). The
mechanistic question of **how** MPs distribute and leach in
microvasculature has only been addressed phenomenologically. To our
knowledge, no published 2-D coarse-grained model combines:

- **Deformable RBCs** (Skalak membrane) and **rigid MPs** in the
  same flow,
- **Fick-type leaching with a finite particle-mass reservoir M_p(t)**,
- **Langmuir surface adsorption** on cell membranes,
- **A realistic stenosis geometry**, and
- **Wall-adhesion kinetics** (catch-slip) for MP-endothelium bonds.

SoftFlow has all five in one engine.

## What the simulation measures

| Quantity | How | Headline relevance |
|---|---|---|
| `dose_bottom_at_stenosis`, `dose_top_at_stenosis` | Sum of scalar `C` in the 4-cell wall band, ±30 cells around the stenosis | Localised endothelial dose |
| `dose_bulk` | `Σ C` outside the wall bands | Comparison baseline |
| `M_p_remaining` per MP | `AdvectionDiffusion.getCapsuleMp(id)` | Per-particle reservoir depletion |
| `M_p_remaining` vs MP **size** | Grouped by `type_label` | Size-dependent leaching |
| MP residence at stenosis | Read positions from VTK frames, integrate dwell time over the stenotic x-window | Particle-trapping signature |
| Cell-free layer (CFL) | Phase-3 `lateral_distribution` (post-process) | Margination of MPs vs RBCs/WBCs |
| Spatial dose map | `AdvectionDiffusion.concentration(0)` final field | Publication figure |

## Methods (for the paper)

- **Lattice Boltzmann method, D2Q9 BGK + Latt-Chopard regularised
  collision** for the fluid (Krüger 2017 §5).
- **Immersed boundary method** with 4-pt Peskin kernel and 2
  iterations of multi-direct forcing (Peskin 2002; Krüger 2012).
- **Skalak constitutive law** for RBC membranes (Skalak et al. 1973).
- **Neo-Hookean** for WBCs; **Hookean rigid** for MPs (`is_rigid=True`).
- **Bell + catch-slip wall adhesion** (Bell 1978; Thomas et al. 2008).
- **Fick leaching** with finite reservoir: `J_leach = k_L (C_eq − C_surf)`,
  `dM_p/dt = −max(J_leach, 0)`.
- **Langmuir adsorption** on RBC/WBC surfaces:
  `dΓ/dt = k_a C (1 − Γ/Γ_max) − k_d Γ`.
- All physics is coarse-grained and lives in 2-D. We are explicit
  about the limitations in the paper's Discussion (see below).

## Dimensionless numbers

| Number | Definition | Value | Regime |
|---|---|---|---|
| Re | `ū H / ν` | ~35 | Inertial-viscous balance for capillary flow |
| Pe | `ū H / D` | ~560 | Advection-dominated chemical transport |
| Da | `k_L L / ū` | ~1.1 | Leaching-diffusion balanced |
| Bi | `k_L R / D` | ~0.6 | Moderate surface resistance |
| Ca | `μ ū / G_s_RBC` | ~0.10 | Significant RBC deformation |

## Limitations to declare

1. **2-D only.** The third dimension matters for true RBC margination
   (~30 % effect, Müller-Fedosov-Gompper 2014). We argue for the
   *trend* and not absolute numbers.
2. **No real chemistry.** "Plasticizer" is a single passive scalar
   with one Fick rate constant. Real BPA / phthalates / PFAS have
   different `k_L`, `D`, and tissue-uptake kinetics.
3. **No endothelial transport.** Wall absorption is treated as a
   sink boundary condition. Penetration of the wall is **out of
   scope**.
4. **Pulsatile flow is approximated as steady.** Project 2 adds
   pulsatility on a related question.
5. **No platelets, no immune response, no clotting.** Pure
   biomechanical transport.

## Running

```bash
python research/01_microplastic_stenosis/run.py    # ~10-15 min
python research/01_microplastic_stenosis/analyse.py # figures
```

Output:

- `vtk_mp_stenosis/{fluid,particles}/*.pvd` — open in ParaView.
- `history.npz` — final scalar field + per-particle quantities.
- `fig1_dose_map.png` — heatmap of concentration with stenosis
  outline (Figure 4 candidate).
- `fig2_per_particle.png` — boxplot of M_p by type (Figure 5
  candidate).
- `fig3_wall_vs_bulk.png` — bar chart of wall vs bulk dose (the
  HEADLINE figure).

## Suggested paper structure

1. **Introduction** — Leslie/Ragusa/Campen MP-in-human findings;
   motivate the stenosis-as-trap hypothesis.
2. **Methods** — LBM-IBM + Skalak RBC + finite-reservoir leaching +
   Langmuir uptake. Brief V&V (Poiseuille profile match, Aris-Taylor
   dispersion).
3. **Results** —
   - Spatial dose map (Fig 1).
   - Wall vs bulk dose ratio for varying stenosis severity (sweep).
   - M_p depletion vs particle size (Fig 2).
   - Margination index for each species (Phase-3 metric).
4. **Discussion** — Limitations, comparison with experimental data
   where possible, and outlook for 3-D extension.
5. **Conclusions** — Stenosis is a *kinetic amplifier* for MP dose
   at the endothelium.

## Target journals

| Journal | IF (2024) | Fit |
|---|---|---|
| *Lab on a Chip* | ~7 | Strong — microfluidic methods + clinical relevance |
| *Environmental Science & Technology* | ~11 | Strong — directly addresses MP-health impact |
| *Soft Matter* | ~3.5 | Solid — coarse-grained mechano-chemistry |
| *Biophysical Journal* | ~3.6 | Possible if methods section is rigorous |
