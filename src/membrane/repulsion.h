#pragma once
#include "../core/parameters.h"
#include "../geometry/obstacle.h"
#include "capsule.h"
#include "capsule_system.h"
#include "cell_list.h"

namespace softflow {

class RepulsionForce {
public:
    explicit RepulsionForce(const RepulsionParams& params) : params_(params) {}

    /// Set domain width for periodic minimum-image in x.
    /// Call with nx > 0 to enable, 0 to disable (non-periodic).
    void setPeriodicX(Real Lx) { Lx_ = Lx; }

    void computeInterCapsule(Capsule& ci, Capsule& cj);
    void computeWallRepulsion(Capsule& c, Real y_bottom, Real y_top);
    void computeAll(CapsuleSystem& system, Real y_bottom, Real y_top);
    void computeObstacleRepulsion(CapsuleSystem& system, const Obstacle& obs);

    /// One-sided repulsion: accumulate force on target from source only.
    /// Thread-safe when each thread owns a unique target capsule.
    void computeOneSidedRepulsion(Capsule& target, const Capsule& source);

private:
    RepulsionParams params_;
    Real Lx_ = 0.0;  ///< domain width (0 = non-periodic)
    CellList cell_list_;  ///< spatial acceleration structure

    /// Apply minimum-image convention in x if periodic
    inline Real minImageDx(Real dx) const {
        if (Lx_ > 0.0) {
            if (dx >  0.5 * Lx_) dx -= Lx_;
            if (dx < -0.5 * Lx_) dx += Lx_;
        }
        return dx;
    }

    /// Cell-list accelerated inter-capsule repulsion (one-sided, thread-safe).
    /// For each membrane node of capsule i, checks only nodes in the same
    /// and neighboring cells rather than all nodes of all other capsules.
    void computeAllCellList(CapsuleSystem& system, Real y_bottom, Real y_top);
};

} // namespace softflow
