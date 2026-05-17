// V&V: Taylor dispersion (CLAUDE.md §6, row 6)
//
// In a 2D plane Poiseuille flow with mean velocity ū driving advection
// of a passive scalar with molecular diffusivity D, a localized pulse
// disperses with an effective long-time diffusivity
//
//   D_eff = D + Pe² D / 210,        Pe = ū H / D,
//
// (Aris 1956 / Taylor 1953). We verify the second moment of the
// concentration distribution along x grows linearly in time with slope
// 2 D_eff after the cross-channel mixing time τ_mix = H² / (4 D).
//
// Tolerance per CLAUDE.md §6: < 5%.
//
// This test is the most demanding regression for the Phase-1 ADR
// rewrite — a first-order ADR equilibrium fails this check by 30+%
// at the Pe used here.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "engine/simulation.h"
#include "lbm/advection_diffusion.h"
#include "lbm/lattice_field.h"

#include <cmath>

using namespace softflow;

TEST_CASE("Taylor dispersion: D_eff = D + Pe^2 D / 210, < 5% error") {
    SimulationParams p;
    // A wide channel keeps Pe meaningful while limiting cost. Tube
    // half-width a = (ny-2)/2 = 11 in lattice units. Driving body force
    // produces ū ~ 0.005, so with D = 0.02 we get Pe = ūH/D ~ 6.
    p.nx = 600;
    p.ny = 24;
    p.dt = 1.0;
    p.fluid.tau           = 0.6;          // ν = 1/30 ≈ 0.0333
    p.fluid.rho0          = 1.0;
    p.fluid.boundary_type = BoundaryType::PERIODIC;
    p.fluid.periodic_y    = false;
    p.fluid.body_force_x  = 5.0e-6;
    p.fluid.collision_model = CollisionModel::BGK;

    p.scalar.enabled              = true;
    p.scalar.n_species            = 1;
    p.scalar.diffusivity          = {0.02};
    p.scalar.inlet_concentration  = {0.0};
    p.scalar.periodic_y           = false;

    p.enable_stability_checks = false;
    p.enable_profiling        = false;
    p.metrics_interval        = 0;
    p.vtk_dump_every = p.csv_dump_every = p.probe_dump_every = p.stats_dump_every = 0;
    p.output_dir = "/tmp/softflow_vv_taylordisp";

    Simulation sim(p);
    sim.initialize();

    // Spin-up the velocity profile to fully developed Poiseuille.
    const int spinup = 30000;
    for (int s = 0; s < spinup; ++s) sim.step();

    // Compute the realised mean velocity by integrating u_x across the
    // fluid cross-section (excluding the two wall sites).
    auto& field = sim.lbmSolver().field();
    Real u_mean = 0.0;
    int  count = 0;
    for (int y = 1; y < p.ny - 1; ++y) {
        for (int x = 0; x < p.nx; ++x) {
            u_mean += field.getUx(x, y);
            ++count;
        }
    }
    u_mean /= count;
    REQUIRE(u_mean > 1e-6);

    const Real D    = p.scalar.diffusivity[0];
    const Real H    = static_cast<Real>(p.ny - 2);
    const Real Pe   = u_mean * H / D;
    const Real D_eff_expected = D + Pe * Pe * D / 210.0;

    // Seed a narrow Gaussian pulse spanning the channel.
    auto* adr = sim.advectionDiffusion();
    REQUIRE(adr != nullptr);
    Real* C = adr->concentrationData(0);
    const Real x0 = p.nx / 4.0;
    const Real sigma0 = 6.0;
    for (int y = 1; y < p.ny - 1; ++y) {
        for (int x = 0; x < p.nx; ++x) {
            Real dx = x - x0;
            C[y * p.nx + x] = std::exp(-dx * dx / (2.0 * sigma0 * sigma0));
        }
    }
    // Re-equilibrate distributions to match.
    {
        Real* g = adr->gData(0);
        static constexpr Real w[9] = {
            4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0,
            1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
        for (int n = 0; n < p.nx * p.ny; ++n) {
            for (int q = 0; q < 9; ++q) g[n * 9 + q] = w[q] * C[n];
        }
    }

    // Compute cross-section-averaged C(x) and its variance about its
    // mass centre. Comoving frame removes ū·t drift.
    auto sectionMoments = [&](Real& mass, Real& mean, Real& var) {
        std::vector<Real> Cbar(p.nx, 0.0);
        for (int x = 0; x < p.nx; ++x) {
            for (int y = 1; y < p.ny - 1; ++y) Cbar[x] += C[y * p.nx + x];
            Cbar[x] /= (p.ny - 2);
        }
        mass = 0.0; mean = 0.0;
        for (int x = 0; x < p.nx; ++x) { mass += Cbar[x]; mean += Cbar[x] * x; }
        mean /= mass;
        var = 0.0;
        for (int x = 0; x < p.nx; ++x) {
            Real dx = x - mean;
            var += Cbar[x] * dx * dx;
        }
        var /= mass;
    };

    // Run past the Aris–Taylor cross-channel mixing time τ_mix = H²/(4D).
    // Then sample variance growth between two instants in the long-time
    // regime.
    const int t_mix    = static_cast<int>(H * H / (4.0 * D));
    const int t_start  = 5 * t_mix;
    const int t_window = 5 * t_mix;

    for (int s = 0; s < t_start; ++s) sim.step();
    Real mass1, mean1, var1; sectionMoments(mass1, mean1, var1);
    for (int s = 0; s < t_window; ++s) sim.step();
    Real mass2, mean2, var2; sectionMoments(mass2, mean2, var2);

    const Real D_eff_meas = (var2 - var1) / (2.0 * t_window);
    const Real rel_err    = std::abs(D_eff_meas - D_eff_expected) / D_eff_expected;

    INFO("ū = " << u_mean << ", Pe = " << Pe);
    INFO("D_eff expected = " << D_eff_expected);
    INFO("D_eff measured = " << D_eff_meas);
    INFO("relative error = " << rel_err);
    CHECK(rel_err < 0.05);
}
