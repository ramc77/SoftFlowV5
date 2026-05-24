# Project 4 — Shear-activated drug carriers in a stenotic microchannel

**A 2D LBM-IBM phase diagram of carrier release-threshold versus
stenosis severity. Engineered shear-activated nanotherapeutics
release their payload preferentially in the elevated-shear zone
inside a vessel constriction; this sweep finds the (γ̇_th, severity)
sweet spot for maximum drug deposition at the stenosis wall.**

## Headline title

> **A 2D LBM-IBM phase diagram of shear-activated drug-carrier
> deposition in stenotic microchannels: tuning release thresholds
> against vessel-narrowing severity.**

## Why this project (and why not the previous ones)

The bidisperse rearrangement attempt (now-deleted Project 4 v1)
showed conclusively that 2D LBM-IBM cannot produce robust
margination — finite-N noise dominated the multi-seed phase
diagram. The drug-penetration-to-tumour attempt (Project 5)
failed because the tumour cluster drifts away from the
released-drug field, leaving an end-of-run snapshot empty.

This project deliberately avoids both failure modes:

| Old failure | How this project avoids it |
|---|---|
| Bidisperse margination → finite-N noise | Single carrier population vs *static* wall obstacle; no multi-particle statistical ensemble |
| Drug-tumour co-location → snapshot-fragile | The metric is **cumulative wall absorption over the run** — monotone-growing, intrinsically time-averaged |
| Tumour-cluster drift confounds delivery | Stenosis is a **rigid obstacle** that does not move; release zone and absorption zone are co-located by construction |
| Mode comparison requires identical equilibration | All cells reach steady-state deposition trajectories; the phase-diagram observable is the *deposition rate*, not an end-of-run snapshot |

## Research question

> **For an engineered shear-activated drug carrier flowing through a
> stenotic microchannel, how should the carrier's release-activation
> threshold $\dot\gamma_\mathrm{th}$ be tuned against the
> stenosis severity (channel-narrowing fraction $\sigma$) to
> maximise local drug deposition while minimising off-target loss?**

The clinical motivation is direct: shear-activated nanotherapeutics
(Korin et al., *Science* **337**, 738, 2012) are an active area of
nanomedicine R&D for treating obstructed blood vessels (stroke,
thrombosis, atherosclerotic plaques). The carriers are designed to
release their payload only above some critical local shear rate —
exactly the rate found at vessel constrictions. The unsolved design
question is the optimal $\dot\gamma_\mathrm{th}$ for a given target
geometry. No systematic 2D LBM-IBM phase diagram of
$(\dot\gamma_\mathrm{th}, \sigma)$ exists in the literature.

## Why publishable now

| Reference | Relevance |
|---|---|
| N. Korin et al., *Science* **337**, 738 (2012) | Foundational shear-activated nanotherapeutic paper |
| S. Mitragotri, "Innovations in drug-delivery devices and biomaterials," *Nat. Rev. Mater.* **9**, 235 (2024) | 2024 review highlighting the design-rule gap |
| P. Decuzzi, M. Ferrari, *Pharm. Res.* **25**, 1399 (2008) | Vascular-carrier design principles |
| A. R. Pries, T. W. Secomb, P. Gaehtgens, *Cardiovasc. Res.* **32**, 654 (1996) | Microvascular hemodynamics |
| W. Higuchi, *J. Pharm. Sci.* **56**, 315 (1967) | Diffusion-controlled release foundation |
| J. Siepmann, F. Siepmann, *J. Control. Release* **161**, 351 (2012) | Release-kinetics taxonomy |
| T. Krüger et al., *The Lattice Boltzmann Method* (Springer, 2017) | LBM-IBM engine basis |

The gap: a clean (release-threshold × stenosis-severity) phase
diagram with a quantitative design rule (i.e. "for severity X, use
threshold Y for maximum deposition").

## What's measured (per cell)

| Field | Definition | Role |
|---|---|---|
| `eta_deposit` | $M_\mathrm{absorbed,target} / M_\mathrm{released}$ at end of run | **Headline targeting efficiency** |
| `eta_offtarget` | $M_\mathrm{absorbed,off}/M_\mathrm{released}$ | Off-target loss diagnostic |
| `eta_remaining_in_fluid` | $1 - \eta_\mathrm{deposit} - \eta_\mathrm{off}$ | Drug still in the fluid at end of run |
| `release_fraction` | $M_\mathrm{released}/M_\mathrm{loaded}$ | Carrier release activity |
| `t_50_deposit` | First step at which $\eta_\mathrm{deposit}(t) \ge \tfrac{1}{2}\eta_\mathrm{deposit,final}$ | Delivery timing |
| `wall_shear_peak` | $\dot\gamma_\mathrm{max}$ at the stenosis throat | Geometry verification |
| `release_fraction_in_throat` | Fraction of released payload that fired inside the throat region | Direct spatial-selectivity measure |

All metrics are **monotone-growing time integrals**, so the noise
mode that killed the bidisperse project (finite-N positional
fluctuation) doesn't apply here.

## Sweep design

**Two-axis grid** — no symmetry, both axes are distinct
physical quantities:

| Stenosis severity $\sigma$ | Throat width (LU) | Geometric reference |
|---|---|---|
| **0 %** | 80 (full) | no-stenosis control; expect zero shear gradient |
| **25 %** | 60 | mild |
| **50 %** | 40 | moderate (physiological-arterial scale) |
| **75 %** | 20 | severe (near-occlusive) |

| Shear-activation threshold $\dot\gamma_\mathrm{th}$ |
|---|
| 5 × 10⁻⁴ (low — sigmoid fires almost everywhere in the channel) |
| 1 × 10⁻³ (medium — fires only when shear is above bulk Poiseuille) |
| 2 × 10⁻³ (high — fires only in the stenosis-induced high-shear band) |
| 5 × 10⁻³ (very high — almost never fires) |

**4 × 4 = 16 cells.** Each cell runs for 30 000 production steps
after a 2 000-step warm-up. Wall-time estimate: ~6 min per cell on
a 16-core laptop, **~1.5 h for the full sweep**.

Smoke version (`--smoke`): 2 × 2 sub-grid × 5 000 steps in ~10 min.

## Design choices

- **Channel geometry**: 400 × 80 LU, periodic in $x$, no-slip walls
  in $y$. Same proven-stable geometry as the DLD project.
- **Stenosis obstacles**: two circular obstacles of radius
  $R_p = 50$ LU centred just outside the top and bottom walls,
  shifted inward to produce the desired throat width.
  $\sigma = 0$ corresponds to no obstacles; $\sigma = 75 \%$
  protrudes each obstacle 30 LU into the channel (throat = 20 LU).
- **Carriers**: 18 stiff Skalak capsules
  ($R = 2$ LU, $G_s = 0.5$), seeded in the upstream region
  $x \in [10, 80]$. Periodic-$x$ → carriers re-circulate past the
  stenosis multiple times during a 30 000-step run.
- **Release kinetic**: `ShearTriggered`, sigmoid
  $k_\mathrm{eff}(\dot\gamma) = k_\mathrm{max}/(1 + \exp(-\beta(\dot\gamma - \dot\gamma_\mathrm{th})))$
  with $k_\mathrm{max} = 0.01$, sigmoid sharpness $\beta = 5 \times 10^3$.
- **Drug field**: single scalar species, diffusivity $D = 0.05$,
  initial concentration 0.
- **Target absorbers**: two narrow `WallAbsorber` bands placed
  **immediately adjacent** to the stenosis obstacles at the throat
  (one above bottom obstacle, one below top obstacle).
  First-order absorption with rate $k_\mathrm{abs} = 0.1$.
- **Off-target absorbers**: two `WallAbsorber` bands placed against
  the channel walls *away* from the stenosis ($x \in [0, 100]$ and
  $x \in [300, 400]$, full $y$-span minus 2 LU near the walls).
  Same first-order rate.
- **Body force**: $F_x = 1.5 \times 10^{-5}$, giving an empty-channel
  Poiseuille centre-line velocity $u_\mathrm{max} \approx 0.014$ LU,
  Mach number $\approx 0.024$, Reynolds $\approx 25$.

## Predicted findings

1. **Diagonal sweet-spot ridge**: maximum $\eta_\mathrm{deposit}$
   along a curve where $\dot\gamma_\mathrm{th}$ matches the
   throat shear rate (which itself grows with severity).
2. **$\eta_\mathrm{deposit}$ collapses** at no-stenosis cells
   (no shear elevation → no release → no deposition).
3. **$\eta_\mathrm{offtarget}$** is highest at low $\dot\gamma_\mathrm{th}$,
   independent of severity (drug released everywhere, captured by
   the upstream/downstream off-target zones).
4. **`t_50_deposit`** is shortest at the sweet-spot cells — the
   delivery happens efficiently per pass.

The phase diagram of $\eta_\mathrm{deposit}$ across the 16 cells
is the headline figure. The sweet-spot ridge is the engineering
design rule.

## Limitations to declare

- **2D only.** Real stenoses are 3D (axisymmetric for arteries,
  asymmetric for atherosclerotic plaques); the 2D simulation
  captures the *geometric narrowing* effect but loses
  cross-stream lateral migration.
- **First-order wall absorption.** A real vessel-wall sink is
  more complex (receptor density, internalisation kinetics).
  We make the simplest assumption for the methodological mapping.
- **Single carrier stiffness.** Carrier deformability under high
  shear (in the throat) may itself modulate release kinetics
  in real carriers; we hold $G_s$ fixed to isolate the threshold
  effect.
- **Steady-state metric.** Drug that exits via off-target absorption
  or remains in fluid at end-of-run is *not* the same as drug that
  *would* eventually deposit — the 30 000-step window is a finite
  snapshot. We report the time-resolved $\eta_\mathrm{deposit}(t)$
  in the per-cell deep-dive.

## Running

```bash
# Smoke (~10 min): 2×2 grid × 5000 steps
python research/04_shear_activated_stenosis/sweep.py --smoke

# Production (~1.5 h): 4×4 grid × 30 000 steps
python research/04_shear_activated_stenosis/sweep.py

# Single-cell baseline (Figure 1 of paper)
python research/04_shear_activated_stenosis/run.py

# Status check (always free)
python research/04_shear_activated_stenosis/status.py

# Figures + summary after sweep
python research/04_shear_activated_stenosis/sweep_analyse.py
```

Resume-safe (same cache pattern from earlier projects):
`Ctrl-C` mid-cell loses only that cell; cached cells skipped on
re-launch. `--force` overrides cache.

## Output layout

```
research/04_shear_activated_stenosis/
├── sweep_out/
│   ├── sev00_th1/                       per-cell directory
│   │   ├── vtk_stenosis/...             ParaView trajectory + scalar field
│   │   ├── history.npz                  per-step deposition + release time series
│   │   ├── cell_result.json             cached summary
│   │   └── config/run_manifest.json     Phase-1 provenance
│   ├── sev00_th2/
│   ├── sev25_th1/
│   ├── ...
│   └── sev75_th4/
├── sweep_results.npz                    aggregated 16-cell grid
└── after `sweep_analyse.py`:
    ├── fig_eta_deposit.png              [HEADLINE] η_deposit × (severity, th)
    ├── fig_eta_offtarget.png            off-target loss heatmap
    ├── fig_release_fraction.png         carrier release activity heatmap
    ├── fig_t50_deposit.png              log10(t_50_deposit) heatmap
    ├── fig_sweet_spot.png               η_deposit vs threshold per severity (story plot)
    └── fig_sweep_summary.txt            tabular summary
```

## Suggested paper structure

1. **Introduction** — Korin 2012 *Science* foundational result;
   the design-rule gap (which threshold for which stenosis);
   prior 3D simulations vs the 2D LBM-IBM gap; engineering value.
2. **Methods** — LBM-IBM (re-use the equation set from the DLD
   paper); stenosis geometry; carrier model; shear-triggered
   release sigmoid; target/off-target absorber definition;
   sweep design.
3. **Results**
   - Single-cell anatomy: snapshots of one representative
     sweet-spot cell showing carrier-deposition cascade (Figure 1).
   - Time series: drug field at three time points + cumulative
     deposition $\eta_\mathrm{deposit}(t)$ (Figure 2).
   - Headline phase diagram: $\eta_\mathrm{deposit}$ on the
     (severity × threshold) grid (Figure 3).
   - The sweet-spot ridge (Figure 4 — story plot, one curve
     per severity).
   - Off-target loss heatmap (Figure 5 — quantifies the cost
     of mis-tuned thresholds).
4. **Discussion** — Mechanism: high-threshold carriers achieve
   spatial selectivity but require severe stenosis to fire;
   trade-off with delivery speed; comparison with Korin 2012
   experimental observation; limitations.
5. **Conclusion** — Design rule:
   $\dot\gamma_\mathrm{th}^\star \approx \dot\gamma_\mathrm{throat}(\sigma)$.
   Future work: 3D extension; multi-species carriers; pulsatile flow.

## Target journals

| Journal | IF (2024) | Fit |
|---|---|---|
| *Microfluidics and Nanofluidics* | ~3.2 | Strong — companion to your DLD paper |
| *Journal of Controlled Release* | ~10.5 | Primary — engineered carrier focus, high impact |
| *Acta Biomaterialia* | ~9.7 | Strong — biomaterial drug delivery |
| *Lab on a Chip* | ~7.0 | If framed as chip-design rule |
| *Phys. Rev. Fluids* | ~3.0 | Backup — capsule mechanics in stenotic flow |
