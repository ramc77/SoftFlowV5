# SoftFlow research pipelines

Three publishable research projects that exploit SoftFlow's unique
capability stack to fill gaps in the 2023-2026 soft-matter-in-flow
literature. Each pipeline is a self-contained directory with a
runnable simulation, an analysis script, a README, and a focused
literature anchor list.

## What makes SoftFlow well-suited

By the end of Phase 5 the codebase combines, in a single 2-D LBM-IBM
engine:

- **Multiple deformable particle types** with selectable membrane
  laws (Hookean / Neo-Hookean / Skalak / WLC), polydisperse insertion,
  and Bell + catch/slip adhesion.
- **Multi-species advection-diffusion** with Fick leaching, finite
  particle-mass reservoirs, and Langmuir surface adsorption.
- **Phase-3 diagnostics** that go beyond bulk margination: Lacey/
  Danckwerts mixing indices, per-type contact matrix Z_ij, lane-order
  parameter, Hoshen-Kopelman clustering, force-network percolation,
  Falk-Langer D²_min.
- **Reproducibility manifest** (git SHA, compiler flags, fully
  resolved params, RNG seed) on every run.

This combination is rare. Most existing codes have either fluid +
rigid particles (lattice Boltzmann), or biological membranes without
chemistry (3-D capsule codes), or chemistry without proper bond
kinetics. We have all three in one engine — that's the publication
leverage.

## The three projects

### Project 1 — `01_microplastic_stenosis/`

> **Question.** Does vessel stenosis preferentially trap microplastics,
> and does this trapping concentrate plasticizer leaching at the
> stenotic wall?

**Topical for 2024-2026.** Microplastics have been detected in human
blood (Leslie et al., *Environ Int* 2022), placenta (Ragusa et al.,
*Environ Int* 2021), and brain tissue (Campen et al., *Tox Sci* 2023).
The transport mechanism in microvasculature — and the spatial
distribution of leached chemicals — has not been simulated with
deformable RBCs + Fick leaching + Langmuir surface adsorption in a
realistic stenosed geometry.

**SoftFlow features used.** Skalak RBCs + Neo-Hookean WBCs + rigid
microplastics; Fick leaching with finite mass reservoir; Langmuir
adsorption on RBC surfaces; circular obstacle pair for the stenosis;
Phase-3 spatial dose map; Phase-1 run_manifest for traceability.

**Target journals.** *Lab on a Chip*, *Environmental Science &
Technology*, *Soft Matter*.

---

### Project 2 — `02_ctc_cluster_pulsatile/`

> **Question.** Does physiological pulsatility (heart-beat-like
> sinusoidal shear) stabilise or fragment circulating-tumour-cell
> (CTC) clusters compared to steady shear, given the catch-slip
> nature of cell-cell adhesion bonds?

**Topical for 2023-2025.** Aceto et al., *Cell* 2014, showed CTC
clusters are the dominant driver of metastasis. Catch-slip bonds
(Thomas-Vogel-Sokurenko, *Annu Rev Biophys* 2008) are well-
characterised at constant force but their behaviour under *time-
varying* force has never been directly probed in a 2D coarse-grained
simulation that also tracks cluster size + spanning + flow-rate-drop
events.

**SoftFlow features used.** Phase-5 division (kept off here —
testing stability, not growth) + Bell/catch-slip adhesion + cluster
identification via Hoshen-Kopelman + embolization detection +
**pulsatile body force** wired through a step callback.

**Target journals.** *Biophys J*, *Phys Rev Fluids*, *PNAS*.

---

### Project 3 — `03_deformability_dld/`

> **Question.** In a deterministic-lateral-displacement (DLD) pillar
> array, does the critical separation size depend on capsule
> stiffness contrast? Can soft particles be sorted by `G_s` alone
> with no size difference?

**Topical for 2020-2025.** Holm et al., *Lab Chip* 2011 derived the
DLD critical-size formula assuming **rigid** spheres. D'Avino &
Maffettone, *J Fluid Mech* 2015 reviewed soft-particle effects.
Recent work (Henon 2017, Vahidkhah & Bagchi 2015) shows softness
shifts the critical bin, but a systematic 2D study mapping
`(G_s_small, G_s_large)` → sorting fidelity is missing.

**SoftFlow features used.** Multi-obstacle pillar array (10×4 grid
of circles offset row-by-row), bidisperse stiff/soft capsules with
**identical size**, lane-order parameter, per-type contact matrix.

**Target journals.** *Lab on a Chip*, *Microfluidics and
Nanofluidics*.

---

## How to run

Each project is fully self-contained. From the repo root:

```bash
# 1. Microplastic margination in stenosed vessel
python research/01_microplastic_stenosis/run.py

# 2. CTC cluster under pulsatile shear
python research/02_ctc_cluster_pulsatile/run.py

# 3. DLD deformability sorting
python research/03_deformability_dld/run.py
```

Each writes ParaView-readable `.pvd` files plus a `.npz` of the
relevant time-series. Each has an `analyse.py` that produces
publication-quality figures.

## Production runs vs the included parameters

Every script ships with `N_STEPS` low enough that a smoke run
finishes in 5–15 minutes on a laptop. For paper-grade statistics,
raise `N_STEPS` to ~500 000 (multi-hour); for the parameter sweeps
referenced in each README, use a cluster or run them overnight.

The reproducibility manifest is identical between smoke and full-
scale runs, so reviewer-grade traceability is preserved.
