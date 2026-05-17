#include "regularized_bgk.h"
#include "lattice.h"
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace softflow {

RegularizedBGK::RegularizedBGK(Real tau)
    : tau_(tau), omega_(1.0 / tau)
{}

void RegularizedBGK::collide(LatticeField& field) {
    const int N = field.size();
    const Real omega = omega_;

    // SoA: gather pointers to each q-plane
    Real* fp[D2Q9::Q];
    for (int q = 0; q < D2Q9::Q; ++q) {
        fp[q] = field.fData() + q * N;
    }

    const Real* __restrict__ rhoPtr = field.rhoData();
    const Real* __restrict__ uxPtr  = field.uxData();
    const Real* __restrict__ uyPtr  = field.uyData();
    const Real* __restrict__ FxPtr  = field.FxData();
    const Real* __restrict__ FyPtr  = field.FyData();
    const CellType* __restrict__ flags = field.cellTypeData();

#ifdef _OPENMP
    #pragma omp parallel for schedule(static)
#endif
    for (int n = 0; n < N; ++n) {
        if (flags[n] != CellType::FLUID) continue;

        Real rho = rhoPtr[n];
        Real ux  = uxPtr[n];
        Real uy  = uyPtr[n];
        Real Fx  = FxPtr[n];
        Real Fy  = FyPtr[n];

        // 1. Compute equilibrium and gather distributions from SoA
        Real feq[D2Q9::Q];
        Real f_local[D2Q9::Q];
        for (int q = 0; q < D2Q9::Q; ++q) {
            feq[q] = D2Q9::feq(q, rho, ux, uy);
            f_local[q] = fp[q][n];
        }

        // 2. Compute non-equilibrium stress tensor
        Real Pi_xx = 0.0, Pi_xy = 0.0, Pi_yy = 0.0;
        for (int q = 0; q < D2Q9::Q; ++q) {
            Real fneq = f_local[q] - feq[q];
            Pi_xx += fneq * D2Q9::cx[q] * D2Q9::cx[q];
            Pi_xy += fneq * D2Q9::cx[q] * D2Q9::cy[q];
            Pi_yy += fneq * D2Q9::cy[q] * D2Q9::cy[q];
        }

        // 3. Regularize and collide
        constexpr Real cs2 = 1.0 / 3.0;
        constexpr Real cs4_2 = 2.0 * cs2 * cs2;

        for (int q = 0; q < D2Q9::Q; ++q) {
            Real Qxx = D2Q9::cx[q] * D2Q9::cx[q] - cs2;
            Real Qxy = D2Q9::cx[q] * D2Q9::cy[q];
            Real Qyy = D2Q9::cy[q] * D2Q9::cy[q] - cs2;

            Real fneq_reg = (D2Q9::w[q] / cs4_2) *
                            (Qxx * Pi_xx + 2.0 * Qxy * Pi_xy + Qyy * Pi_yy);

            // Guo forcing term
            Real eu = D2Q9::cx[q] * ux + D2Q9::cy[q] * uy;
            Real Fi = (1.0 - 0.5 * omega) * D2Q9::w[q] * (
                (D2Q9::cx[q] - ux + eu * D2Q9::cx[q] / cs2) * Fx / cs2 +
                (D2Q9::cy[q] - uy + eu * D2Q9::cy[q] / cs2) * Fy / cs2
            );

            // Scatter back to SoA plane
            fp[q][n] = feq[q] + (1.0 - omega) * fneq_reg + Fi;
        }
    }
}

} // namespace softflow
