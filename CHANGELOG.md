# Changelog

All notable changes to SoftFlow are recorded here. The format is loosely
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and
this project adheres to [Semantic Versioning](https://semver.org/) once
it reaches 1.0.

## [Unreleased] — Phase 5: Tumour-growth / aggregate-formation extensions

CLAUDE.md §7.4 verbatim. Stochastic division kinetics + volume-
exclusion-aware daughter placement + embolization detection (flow-
rate drop time series + spanning-cluster span time series + events),
all sitting on top of the existing C++ Bell + catch/slip adhesion
machinery (Phase 0) and the Phase-3 `hoshen_kopelman` /
`force_percolation` diagnostics.

The C++ engine is **untouched**. Every Phase-5 component is pure
Python, so Phase-5 changes never risk regressing the V&V suite.

### Strong-language constraint (CLAUDE.md §7.4, mandatory)

This is a **coarse-grained mechano-chemical proxy for circulating-
cell aggregation under flow, NOT a validated cancer model**. The
caveat appears verbatim in [`docs/tumor_growth.md`](docs/tumor_growth.md),
in [`examples/05_tumor_growth/README.md`](examples/05_tumor_growth/README.md),
in every module docstring, and in every public-class docstring.

Specifically:
- 2-D coarse-grained capsules, not realistic tumour cells.
- "Tumour" terminology kept for CLAUDE.md alignment; internal
  language uses "aggregate" / "cluster" / "circulating-cell".
- Division is Poisson-clocked stochastic insertion, not any specific
  cell-cycle / mitogen / oncogenic signalling.
- "Stress" is the magnitude of the local rate-of-strain tensor; NOT
  membrane stress, molecular bond, or calibrated mechano-transducer.
- "Nutrient" is a generic scalar species — we don't simulate
  glucose / oxygen / pH chemistry.
- Embolization detection reports *events*, not clinical outcomes.
- No 3-D, no patient-specific geometry (CLAUDE.md §10).

### Added

- **`pysoftflow.tumor_growth` package** under
  [`python/pysoftflow/tumor_growth/`](python/pysoftflow/tumor_growth/):

  - `division.py` — `DivisionKinetic` ABC + `StressNutrientDivision`
    (Poisson rate, γ̇ threshold, nutrient threshold, optional
    lineage cap). Closed-form Poisson firing probability
    `1 − exp(−k_div · dt)`, exact for any dt.

  - `daughters.py` — `DaughterPlacer.propose(...)` tries random
    angles around a ring at `factor · (parent_r + daughter_r +
    min_gap)`. Each candidate is checked against walls, obstacles
    (signed-distance callable), and existing capsules under
    periodic-x minimum-image.

  - `embolization.py` — `EmbolizationDetector` with single-x_section
    flow-rate measurement Q(t) = Σ u_x(x_section, j), Phase-3
    `force_percolation` for spanning-cluster detection, and
    `EmbolizationEvent` records that fire only when *both* gates
    fire simultaneously. Per-step `flow_rate_history` and
    `cluster_span_history` accessors.

  - `runner.py` — `TumorGrowthRun` orchestrator with single per-
    step callback. Snapshot eligible capsule indices at start of
    step (so daughters added now don't divide in the same step),
    probe shear + nutrient at each parent's centroid, call kinetic
    + placer, append daughter via `CapsuleSystem.addCapsule(...)`,
    run all detectors, append a `RunStepRecord`.

- **Showcase** — [`examples/05_tumor_growth/`](examples/05_tumor_growth/):
  - `run.py` — 200×40 narrow channel, body-force-driven, 5 adhesive
    seed capsules with Bell adhesion + `--catch-slip` toggle, 1.0
    uniform nutrient prefill, embolization watch at x=140. Outputs
    VTK + `history.npz` + `events.npz`.
  - `analyse.py` — three matplotlib plots (division history,
    cluster size, flow rate / cluster span overlay).
  - `README.md` — strong-language Limitations block (mandatory).

- **38 new Python tests** across 5 files:
  - `test_tumor_division.py` (10) — threshold gates, 4σ Poisson rate
    verification, lineage cap, seeded reproducibility.
  - `test_tumor_daughters.py` (10) — free-space, narrow-channel
    rejection, dense-neighbourhood rejection, periodic-x wrap,
    obstacle blocking, scaled daughter, reproducibility.
  - `test_tumor_embolization.py` (9) — full gate combinations,
    history accumulation, single-x flow integration, validation.
  - `test_tumor_runner.py` (6) — end-to-end against a live
    `Simulation`, high-stress block, low-nutrient block, monotonic
    division total, idempotent attach, embolization-detector
    pipeline.
  - `test_tumor_example.py` (3) — example default, catch/slip,
    analyse-script round-trip.

- **`docs/tumor_growth.md`** — equations, references (Bell 1978,
  Thomas-Vogel-Sokurenko 2008, Aceto 2014, Au 2016, Stauffer-
  Aharony 1994, Hoshen-Kopelman 1976), and the strong-language
  constraint repeated in full.

### Real find along the way

The default `DaughterPlacer.ring_radius_factor=1.05` initially
applied to `(parent_r + daughter_r)`, which placed daughters at
distance `1.05·(r_p + r_d) = 4.2` from the parent — *inside* the
parent's `r_p + r_d + min_gap = 4.5` exclusion envelope. The placer
then rejected every candidate, regardless of angle. Fixed by
including `min_gap` in the ring radius: `factor · (r_p + r_d +
min_gap)`. The unit tests caught it on the very first run; the
showcase example would have produced zero divisions silently
otherwise.

### Known follow-ups (in `docs/tumor_growth.md`)

- Per-capsule strain accessor (membrane stress proxy instead of
  fluid shear). Small C++ binding addition.
- Daughter inherits zero bonds; partial inheritance is a follow-up.
- Single-x flow rate vs. multi-section averaged Q(t).
- Daughter "growing" from r → r_target over a few steps for softer
  initial physics.
- 3-D, patient-specific geometry (CLAUDE.md §10 out-of-scope).

---

## [Unreleased] — Phase 4: Drug-delivery extensions

CLAUDE.md §7.3 verbatim. Five release-kinetic models for carriers,
wall-region absorbers (first-order and Michaelis-Menten), four
delivery metrics, and a stiffness × release-rate sweep showcase.

The C++ engine is **untouched**. Every Phase-4 component is pure
Python on top of the existing chemistry path (Phase-1's Fick leaching
+ finite reservoir M_p + Langmuir adsorption + Peskin spreading), so
Phase-4 changes never risk regressing the V&V suite.

### Added

- **`pysoftflow.drug_delivery` package** under
  [`python/pysoftflow/drug_delivery/`](python/pysoftflow/drug_delivery/):

  - `kinetics.py` — five release models with a uniform
    `update(carrier_state, fluid_probe, dt) -> released_mass` API:
    - `DiffusionControlled(k_leach, C_eq)` — adapter to the existing
      C++ Fick path (Higuchi 1961).
    - `FirstOrder(k_rel)` — exact closed-form decay
      `ΔM = M (1 − exp(−k_rel · dt))`.
    - `ShearTriggered(k_max, γ_thresh, sharpness)` — sigmoid in the
      magnitude of the local rate-of-strain tensor at the carrier
      centroid (Bao & Suresh 2003 in spirit).
    - `PhTriggered(k_max, C_thresh, species, sharpness)` — sigmoid
      in a designated trigger species. Honest naming: not pH
      chemistry, just second-scalar-triggered (Schmaljohann 2006).
    - `Burst(release_time, fraction)` — one-shot drop, per-carrier
      independent firing.

  - `absorbers.py` — `WallAbsorber` rectangular patch with first-
    order or Michaelis-Menten kinetics. Operates on the writable
    numpy view of the C++ scalar field; tracks
    `cumulative_absorbed` + per-step `history` array.

  - `metrics.py` — `delivery_efficiency`, `off_target_fraction`,
    `residence_time_distribution` (with type filter and histogram),
    `spatial_dose_map` (with overlap summing).

  - `runner.py` — `DrugDeliveryRun` orchestrator with single
    per-step callback that refreshes M_p from C++, builds a
    `FluidProbe` (concentration + shear rate + species probes),
    drives every kinetic, Peskin-spreads released mass for non-
    Fickian models (3×3 uniform spreader; coarser than the full
    Peskin 4-pt kernel — fine for relative comparisons, see
    Limitations), applies all absorbers, and records a per-step
    `RunStepRecord`.

- **Showcase** — [`examples/04_drug_delivery/`](examples/04_drug_delivery/):
  - `run.py` — single-configuration showcase. CLI `--mode {diffusion,
    first_order, shear, ph, burst}` selects the kinetic. Outputs
    `history.npz` + standard SoftFlow VTK + `run_manifest.json`.
  - `sweep.py` — runs `run.py`'s `build_simulation` over a
    (G_s, k_release) grid; produces `results.npz` + a matplotlib
    heatmap. `--smoke` (2×2), default (3×3), `--full` (5×5).
  - `README.md` — limitations, dimensionless-number table
    (Re, Pe, Da_I, Da_II, Ma, χ), reproducibility instructions.

- **53 new Python tests** (six files):
  - `test_drug_kinetics.py` (17) — closed-form decay, sigmoid
    asymptotics, threshold sharpness, per-carrier independence.
  - `test_drug_absorbers.py` (13) — analytic decay, mass conservation,
    MM saturation/linearity, over-depletion clamp.
  - `test_drug_metrics.py` (14) — η/OTF idealised mass balance,
    RTD with type filter, dose-map overlap behaviour.
  - `test_drug_runner.py` (3) — end-to-end pipeline against a live
    `Simulation`.
  - `test_drug_example.py` (6) — every release mode + sweep
    mass-balance, parametrised CI.

- **`docs/drug_delivery.md`** — algorithmic notes, references
  (Higuchi 1961, Schmaljohann 2006, Bao & Suresh 2003, Michaelis-
  Menten 1913, Decuzzi & Ferrari 2008), and the strong-language
  constraint repeated in full.

### Strong language constraint (per CLAUDE.md §7.4 spirit / §7.3 letter)

This entire module is a **methodological-exploration toolkit**, not
a validated drug-carrier or pharmacokinetic model. The caveats appear
verbatim in `docs/drug_delivery.md`, in `examples/04_drug_delivery/
README.md`, and in every kinetic class's docstring:

- 2-D coarse-grained capsules, not realistic carriers.
- Wall absorbers are boundary-condition sinks, not tissue models.
- "pH-triggered" is sigmoid in a designated species, not real pH.
- "Shear-triggered" is sigmoid in |D|, not any specific bond.
- Non-Fickian spreading is 3×3 uniform, not the full Peskin kernel.

### Known follow-ups

- **Per-Peskin-cell spreading** for non-Fickian release modes (current
  3×3 uniform is fine for relative comparisons, not for quantitative
  cell-by-cell deposition profiles).
- **MM uptake with explicit transporter saturation** (current k_cat
  saturation is the simplest M-M form; a two-state
  catalysis-+-binding model would be more realistic).
- **Per-carrier kinetic state for stimuli-responsive models** —
  ShearTriggered / PhTriggered currently re-evaluate the sigmoid
  every step regardless of carrier history. Adding a hysteresis or
  one-way switch is a small follow-up.
- **3-D extension and patient-specific geometry** are CLAUDE.md §10
  out-of-scope.

---

## [Unreleased] — Phase 3: Segregation, mixing, jamming, pattern-formation diagnostics

CLAUDE.md §7.2 verbatim. Twelve post-processing diagnostics consumed
from saved snapshots / live ``Simulation`` objects, plus a single-file
HDF5 export with the run manifest copied in for traceability.

The C++ engine is **untouched**. All diagnostics are pure Python
(numpy + scipy.spatial / scipy.sparse), so Phase-3 changes never
risk regressing the V&V suite.

### Added

- **`pysoftflow.analysis` package** under `python/pysoftflow/analysis/`,
  organised by family:
  - `snapshot.py` — `SimulationSnapshot` data layer with
    `from_simulation(sim)` (live extraction including adhesion bonds
    via `AdhesionModel.getBonds()`) and `from_arrays(…)` (synthetic).
  - `mixing.py` — `lacey_index`, `danckwerts_intensity`,
    `contact_asymmetry` (Lacey 1954, Danckwerts 1952). Returns a
    `MixingIndex` dataclass with σ², σ²_R, σ²_M components for
    transparency.
  - `rdf.py` — species-resolved `radial_distribution(snap, type_a,
    type_b, …)` using `scipy.spatial.cKDTree` + ±Lx image padding
    for periodic-x wrap.
  - `patterns.py` — `lane_order(snap, axis)` ⟨cos 2θ⟩ on velocities;
    `hoshen_kopelman(snap, …)` connected-component labelling on the
    contact / bond graph (uses `scipy.sparse.csgraph`); and
    `cluster_persistence(labels_a, labels_b)` lag-1 Jaccard
    similarity over same-cluster pairs.
  - `jamming.py` — `packing_field(snap, n_x, n_y)` coarse-grained
    φ(x, y); `contact_number(snap, …)` with `Z_per_particle`;
    `per_type_contact_stats(snap, …)` returning the (T, T) matrix
    `Z_matrix[i, j]` of per-type average contacts (verbatim per the
    user request — "smaller particles average … contacts with their
    own type and … with other type"); `force_percolation(snap, …)`
    cross-channel spanning detection; `mean_squared_displacement(
    snaps, …)` with late-half log-log plateau detection;
    `non_affine_d2min(snap_t0, snap_t1, …)` Falk–Langer 1998 D²_min
    via least-squares affine fit over local neighbour patches.
  - `hdf5_export.py` — `write_diagnostics_h5`, `load_diagnostics_h5`,
    `load_run_manifest` for a single-file `/manifest`, `/snapshots`,
    `/diagnostics` HDF5 layout. Recursive dataclass / nested-dict
    traversal so per-snapshot diagnostic dicts land as nested groups
    (not `repr` strings).

- **`examples/02_bidisperse_segregation/analyse.py`** — runs the
  bidisperse simulation, samples snapshots, computes every Phase-3
  diagnostic on each, runs the time-series MSD / persistence /
  D²_min over the whole snapshot list, and writes a single HDF5
  file with the run's `run_manifest.json` included. `--smoke` mode
  for CI.

- **73 Python tests** (`tests/python/test_analysis_*.py`,
  `test_example_analyse.py`) covering each diagnostic with synthetic
  inputs of known answer:
  - random-uniform mix → Lacey M ≈ 1, Danckwerts I_S ≈ 0
  - stripe-segregated → M ≈ 0, I_S ≈ 1 (along the stripe axis)
  - hex lattice → first RDF peak at lattice spacing, Z̄ → 6 interior
  - aligned velocities → lane-order Φ = 1; perpendicular → −1
  - per-type Z_ij identity check `n_i · Z_ij == n_j · Z_ji`
  - hand-crafted spanning chain → percolation True
  - synthetic random walk → MSD slope ~1; trapped → plateau detected
  - pure shear → Falk–Langer D²_min ≡ 0; random shuffle → > 0
  - HDF5 round-trip preserves snapshots, diagnostics, and manifest

- **`docs/analysis/README.md`** — single comprehensive reference for
  the twelve diagnostics, equations, and references (Lacey 1954,
  Danckwerts 1952, Hoshen-Kopelman 1976, Vissers 2011, O'Hern 2003,
  Falk-Langer 1998, Stauffer-Aharony 1994).

### Real find along the way

- The HDF5 writer initially serialised nested dataclasses (e.g. a
  `MixingIndex` field of a per-snapshot dict) as `repr` strings,
  which silently lost the σ²-components and prevented the analyse
  smoke test from inspecting the diagnostic tree. Fixed by detecting
  dataclass / Mapping values inside `_store_field` and recursing
  into a subgroup. Verified by the round-trip test exercising both
  scalar diagnostics and nested-dataclass diagnostics in the same
  file.

### Known follow-ups

- **MSD periodic unwrap** is per-pair minimum-image. For runs that
  wrap many times, build snapshots from absolute-coordinate CSV
  instead of live `Simulation` objects.
- **D²_min uses a fixed neighbour cutoff** at t₀. Adaptive
  neighbour-set choice for large displacements is a follow-up.
- **Plot / notebook layer** is intentionally absent — diagnostics
  emit numpy arrays and HDF5; users plot in their own notebooks.
- **Auto-collect during a run** (vs post-process from snapshots) is
  not wired. The post-process pattern is easier to debug and re-run;
  auto-collect can land if a user asks.
- **Phase-3 example for *INLET_OUTLET* BCs** — current example uses
  periodic, which matches CLAUDE.md §2 but doesn't exercise the
  conveyor-style dynamics.

---

## [Unreleased] — Phase 2: Particle insertion module

CLAUDE.md §7.1 verbatim: static layouts (square, hex, RSA, Poisson-
disk), polydisperse mixtures, dynamic insertion (constant-flux,
Poisson-stochastic, conveyor), region masks (rect, circle, polygon,
image), and a bidisperse-segregation showcase.

Substrate for Phases 3 (segregation diagnostics), 4 (drug delivery),
and 5 (tumour aggregation).

### Added

- **Insertion module** under `src/core/insertion/`, organised around
  three orthogonal abstract bases:
  - `IRegion` — `RectRegion`, `CircleRegion`, `PolygonRegion`,
    `ImageMaskRegion` (PGM P2/P5 reader).
  - `ISizeDistribution` — `Monodisperse`, `Bidisperse`, `Lognormal`
    (hard-truncated), `UserDiscrete` (radii + weights via
    `std::discrete_distribution`).
  - `IInserter` — `SquareLatticeInserter`,
    `HexagonalLatticeInserter`, `RSAInserter`,
    `PoissonDiskInserter` (Bridson 2007).
  - `IDynamicInserter` — `PoissonStochasticInserter`,
    `ConstantFluxInserter`, `ConveyorInserter`.
  - Shared overlap helpers (`isPlacementValid`, `distance`,
    `minImageDx`) that respect periodic-x, walls, obstacles, and
    `min_gap` consistently.

- **`Simulation::insertCapsules(inserter, mparams, type, …)`** — one-
  shot static fill. Builds the `InsertionContext` from the current
  channel/obstacle/capsule state, drives the inserter, and adds each
  returned `Placement` as a capsule. Returns the count placed.

- **`Simulation::registerDynamicInserter(inserter, mparams, type,
  …)`** — register a per-step dynamic inserter. The simulation
  drains it at step 10b (between the y-clamp and ML/metrics) so new
  capsules see full physics on the next step. Each registration owns
  a deterministic `std::mt19937_64` sub-stream derived from
  `params.rng_seed` via FNV-1a hash of a user-supplied `seed_tag`.

- **`pysoftflow.insertion`** Python facade matching the CLAUDE.md
  §7.1 keyword-argument syntax:
  ```python
  inserter = Inserter.hex_lattice(
      region  = RectRegion(x=(0, 200), y=(10, 90)),
      spacing = 8.0,
      sizes   = SizeDistribution.bidisperse(r_small=2.0, r_large=4.0,
                                             fraction_small=0.5),
      jitter  = 0.1)
  sim.insertCapsules(inserter, rbc_params, type=0)
  ```
  All three abstract bases are bound under shared_ptr holders so
  Python users can compose regions, sizes, and inserters freely.

- **`examples/02_bidisperse_segregation/`** — 50/50 small/large
  bidisperse suspension, hex-lattice fill for the small species and
  RSA fill for the large, periodic body-force-driven channel,
  prints the (Re, Ca, χ_s, χ_L, Ma) tuple at startup. `--smoke`
  flag for ~5-second CI verification.

- **`docs/theory/insertion.md`** — algorithmic notes, references
  (Widom, Talbot–Tarjus–Schaaf, Bridson), and the pieces deliberately
  left for future phases (Python user-PDF callback, area-weighted
  log-normal, INLET_OUTLET conveyor).

- **Tests** — 7 C++ unit tests (`test_insertion_scaffold`,
  `_lattice`, `_rsa`, `_size`, `_region`, `_dynamic`, `_image_mask`)
  covering count, no-overlap, reproducibility, jamming density,
  Bridson invariant, statistical fraction recovery, region
  membership, dynamic rate match, and the end-to-end `Simulation`
  pipeline. Plus 7 Python tests exercising the full facade ↔ C++ ↔
  Simulation round-trip and a smoke test for the showcase example.

### Changed

- **`SimulationParams` Python binding** now exposes `rng_seed` and
  `max_lattice_force` (Phase-1 fields that had no pybind11 export).
  Without this, the Phase-2 reproducibility plumbing was unreachable
  from Python.

### Real find along the way

- **`-ffast-math` interacts with `std::numeric_limits<double>::
  infinity()`** in a hostile way: under AppleClang's
  `-ffinite-math-only` (implied by `-ffast-math`), the optimizer
  collapses an `infinity()` assignment to `0`. The Phase-2
  `PoissonDiskInserter::lastMinSeparation()` initially used
  `infinity()` as a "no pair seen yet" sentinel and silently
  returned 0. Fixed locally by switching to
  `std::numeric_limits<Real>::max()`. The broader problem (the same
  flag also breaks `std::isnan`/`std::isinf` calls in
  `Simulation::checkStability`, `immersed_boundary.cpp`, and
  `cell_list.h`) was previously flagged in REVIEW.md §6.5 and
  remains a Phase-1 follow-up — a `-DSOFTFLOW_DETERMINISTIC=ON` build
  mode that drops `-ffast-math` for V&V runs.

### Known follow-ups

- Python user-PDF callback for `ISizeDistribution`.
- Area-weighted (vs number-weighted) sampling for `Lognormal`.
- `ConveyorInserter` semantics for `INLET_OUTLET` BCs (currently
  count-in-region; should optionally measure total-domain count).
- ImageMaskRegion: 16-bit graymaps and PBM (P1/P4) variants.

---

## [Unreleased] — Phase 1: Code review & hardening

This phase addresses the critical and high-severity findings from
[REVIEW.md](REVIEW.md) §6, which catalogued the gaps that needed to
close before SoftFlow could support a peer-reviewed software paper.

### Added

- **`tests/` directory + doctest scaffolding.** The previous build was
  broken on the default `BUILD_TESTS=ON` because `CMakeLists.txt:49`
  called `add_subdirectory(tests)` on a non-existent directory. The
  build now configures cleanly and ships with five passing tests:
  - `tests/unit/test_lattice` — D2Q9 weight, opposite-direction, and
    equilibrium-moment identities.
  - `tests/unit/test_run_manifest` — manifest-writer regression: drives
    a 5-step LBM-only sim and verifies the schema.
  - `tests/verification/test_poiseuille` — parabolic profile,
    L2 error < 1 % (CLAUDE.md §6 row 1).
  - `tests/verification/test_taylor_green` — exponential decay rate
    < 2 % (row 2).
  - `tests/verification/test_pure_diffusion` — Gaussian variance growth
    σ² = σ₀² + 2 D t < 2 % (row 5).

  `tests/verification/TODO.md` documents the remaining V&V cases:
  Taylor dispersion (parked due to a test-side periodic-wrap bug, not
  a physics regression), Skotheim–Secomb tank-treading, and the
  Segré–Silberberg single-disk equilibrium.

- **`docs/theory/`.** Six concise notes (≤ 2 pages each) tying source
  lines to equations and to CLAUDE.md §8 literature anchors:
  - [lbm.md](docs/theory/lbm.md) — D2Q9 with BGK / MRT / regularized.
  - [ibm.md](docs/theory/ibm.md) — Peskin 4-pt direct + multi-direct.
  - [adr.md](docs/theory/adr.md) — D2Q9 ADR + leaching/Langmuir
    chemistry.
  - [membrane.md](docs/theory/membrane.md) — Hookean / Neo-Hookean /
    Skalak / WLC + bending + area / perimeter / Kelvin-Voigt.
  - [adhesion.md](docs/theory/adhesion.md) — Bell + catch–slip kinetics.
  - [lubrication.md](docs/theory/lubrication.md) — Brenner-style 2-D
    regularized closure.

- **`run_manifest.json`** — replaces the old 17-field
  `simulation_config.json`. Each run now writes a comprehensive
  reproducibility record:
  - Build provenance (git SHA, dirty flag, branch, compiler ID +
    version, fully resolved CXX flags including `-ffast-math`,
    OpenMP thread count, system name and processor).
  - ISO-8601 wall-clock timestamp.
  - Canonical RNG seed (new `SimulationParams::rng_seed`).
  - Configurable IBM-LBM force cap (new
    `SimulationParams::max_lattice_force`, default `0.01`).
  - The complete resolved `SimulationParams` tree, with every enum
    serialized as a string name (so reordering enum values does not
    silently invalidate older runs).

  `tests/unit/test_run_manifest` pins this schema.

- **`pysoftflow.units`** — SI ↔ lattice unit converter and dimensionless-
  number helpers. `LatticeUnits.from_channel_flow(...)` derives
  `dx_si`, `dt_si` from channel geometry, viscosity, and a target τ;
  `velocity_to_lb` / `velocity_to_si` (and friends) round-trip exactly;
  `reynolds`, `capillary`, `peclet`, `confinement` compute dimensionless
  numbers from SI inputs. Mach-number guard (`u_lb < 0.1`) refuses
  unstable setups at construction. Eight pytest round-trip tests in
  `tests/python/test_units.py`.

- **Build / packaging.**
  - Root `LICENSE` (MIT, with PI attribution).
  - Root `pyproject.toml` (scikit-build-core, ruff, mypy, pytest).
  - Phase-1 working agreement: ruff `line-length = 100`,
    mypy `strict`, pytest in `tests/python/`.

- **Force-cap visibility.** `Simulation::checkStability()` now reports
  the per-step count of fluid nodes whose IBM/body-force magnitude
  exceeded `params.max_lattice_force` and was rescaled. Previously
  the rescaling was completely silent.

### Changed

- **C++ standard: 17 → 20.** Aligns the build with CLAUDE.md §4.
  `CMAKE_CXX_EXTENSIONS = OFF` for portability. The full codebase
  compiles clean under AppleClang 21.

- **CompilerFlags.cmake.** All compile options now go through
  `softflow_add_flag()`, which simultaneously calls
  `add_compile_options()` and accumulates the flag list into
  `SOFTFLOW_CXX_FLAGS` for the manifest. Users can audit exactly what
  optimization flags compiled their results.

- **ADR collide/stream split** (`src/lbm/advection_diffusion.cpp`).
  The pre-Phase-1 implementation fused collide and stream into a
  single OpenMP-parallel pass and applied solid bounce-back by
  writing into the live `g_[species]` buffer **inside** the parallel
  region — a real race condition that produced nondeterministic
  results at `OMP_NUM_THREADS > 1`. The replacement is a clean
  three-pass pipeline (`collide → streamWithBC → computeConcentration`)
  with disjoint readers and writers in each pass and a textbook
  halfway bounce-back at solid links.

  The equilibrium remains **first-order in u**
  (`g_eq = w_q C (1 + (e_q · u) / cs²)`, Krüger 2017 §8.5.2): an
  earlier draft of this commit "upgraded" it to second order, which
  fails the Aris–Taylor regression by ~10× because the higher-order
  velocity terms pollute the recovered diffusivity. REVIEW.md §6.1
  has been corrected accordingly.

- **`AdhesionParams::adhesion_matrix` empty default.** Previously,
  an empty matrix meant "all type pairs may bond" — a silent foot-gun
  that produced unintended adhesion when users enabled the model
  without configuring type pairs. The new contract is **empty matrix
  ⇒ no adhesion**. The change is announced loudly the first time the
  model is constructed so old scripts surface the regression
  immediately. Migration: set `params.adhesion.adhesion_matrix[i][j]
  = true` for each pair `(i, j)` of capsule types that should be
  allowed to bond.

- **Interpolated bounce-back periodic-x wrap.** When an obstacle sits
  on the streamwise periodic seam, the BFL `q < ½` upstream-link
  fallback used to silently degrade to halfway bounce-back. It now
  wraps through `params.nx`. `LBMSolver` propagates the periodic-x
  flag to `InterpolatedBounceBack` automatically.

- **Force regularization** (`simulation.cpp:319–336`).
  `Fxy_max = 0.01` is no longer hard-coded; it is now
  `params.max_lattice_force` and is captured in the manifest.

- **`examples/README.md`** rewritten with per-example dimensionless-
  number tables, limitation notes, and reproducibility instructions.

### Fixed

- **OpenMP race in ADR bounce-back** at solid walls — see "Changed"
  above. Phase-1's principal correctness fix.

- **Build break on default settings.** `BUILD_TESTS=ON` (the default)
  no longer fails configure. The `tests` subdirectory is also gated
  on the existence of its own `CMakeLists.txt` so an accidental
  deletion produces a clear warning rather than a configure error.

### Reverted

- **ADR equilibrium order.** A draft of this phase moved the ADR
  equilibrium to second-order in **u** based on a misreading of
  REVIEW.md §6.1. This was reverted after the Taylor-dispersion
  regression test exposed the resulting 10× error in the recovered
  diffusivity. The textbook first-order form was already correct.

### Known follow-ups (deferred to a Phase-1 patch or Phase 2)

- `tests/verification/test_taylor_dispersion.cpp`: re-enable after
  rewriting the test with a comoving frame or `INLET_OUTLET` BCs to
  remove periodic-wrap pollution of the variance estimate.
- `tests/verification/test_capsule_shear.cpp`: Skotheim–Secomb tank-
  treading (CLAUDE.md §6 row 3).
- `tests/verification/test_segre_silberberg.cpp`: rigid-disk radial
  equilibrium (row 4).
- IBM δ-stencil clipping diagnostic at y-walls (REVIEW.md §6.3).
- ADR storage layout: nested `vector<vector<Real>>` →
  flat `AlignedArray` (REVIEW.md §6.4 perf item).
- Default to interpolated bounce-back when curved obstacles are
  registered (REVIEW.md §6.1).
- `-DSOFTFLOW_DETERMINISTIC=ON` build mode that drops `-ffast-math`
  for V&V runs (REVIEW.md §6.5).
- Anti-bounce-back ADR inlet (`docs/theory/adr.md` TODO).
- A `velocity-Verlet` capsule integrator for stiff Skalak coefficients
  (`docs/theory/membrane.md`).
