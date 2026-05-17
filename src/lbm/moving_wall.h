#pragma once
#include "lattice_field.h"
#include "../core/parameters.h"

namespace softflow {

/// Moving wall boundary condition.
/// Modified bounce-back with non-zero wall velocity (Ladd 1994).
/// f_opp(x_fluid, opp[q]) = f(x_solid, q) - 2*w_q*rho*(e_q . u_wall)/cs^2
class MovingWall {
public:
    MovingWall(Real top_wall_velocity = 0.0, Real bottom_wall_velocity = 0.0);

    /// Apply moving wall BC after streaming
    void apply(LatticeField& field);

private:
    Real top_vel_;
    Real bot_vel_;
};

} // namespace softflow
