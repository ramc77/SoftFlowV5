#include "repulsion.h"
#include <cmath>

namespace softflow {

void RepulsionForce::computeInterCapsule(Capsule& ci, Capsule& cj) {
    auto& pos_i = ci.positions();
    auto& pos_j = cj.positions();
    auto& f_i = ci.forces();
    auto& f_j = cj.forces();

    Real r_cut = params_.r_cut;
    Real epsilon = params_.epsilon;
    Real sigma = params_.sigma;
    int power = params_.power;

    for (int a = 0; a < ci.numNodes(); ++a) {
        for (int b = 0; b < cj.numNodes(); ++b) {
            Vec2d diff = pos_i[a] - pos_j[b]; // points from j to i
            diff.x = minImageDx(diff.x);       // periodic image in x
            Real r = diff.norm();

            if (r < 1e-15 || r >= r_cut) continue;

            // Morse-like repulsion: F = epsilon * (sigma/r)^power
            Real ratio = sigma / r;
            Real ratio_pow = std::pow(ratio, power);
            Real F_mag = epsilon * ratio_pow;

            Vec2d r_hat = diff / r; // unit vector from j toward i
            Vec2d F = r_hat * F_mag;

            f_i[a] += F;
            f_j[b] -= F;
        }
    }
}

void RepulsionForce::computeWallRepulsion(Capsule& c, Real y_bottom, Real y_top) {
    auto& pos = c.positions();
    auto& f = c.forces();

    Real r_cut = params_.r_cut;
    Real epsilon = params_.epsilon;
    Real sigma = params_.sigma;
    int power = params_.power;

    for (int i = 0; i < c.numNodes(); ++i) {
        // Bottom wall repulsion (force in +y direction)
        {
            Real r = pos[i].y - y_bottom;
            if (r > 0.0 && r < r_cut) {
                Real ratio = sigma / r;
                Real ratio_pow = std::pow(ratio, power);
                Real F_mag = epsilon * ratio_pow;
                f[i].y += F_mag;
            }
        }

        // Top wall repulsion (force in -y direction)
        {
            Real r = y_top - pos[i].y;
            if (r > 0.0 && r < r_cut) {
                Real ratio = sigma / r;
                Real ratio_pow = std::pow(ratio, power);
                Real F_mag = epsilon * ratio_pow;
                f[i].y -= F_mag;
            }
        }
    }
}

/// One-sided repulsion: force on target from source only.
/// Thread-safe when each thread owns a unique target capsule.
void RepulsionForce::computeOneSidedRepulsion(Capsule& target, const Capsule& source) {
    auto& pos_t = target.positions();
    const auto& pos_s = source.positions();
    auto& f_t = target.forces();

    Real r_cut = params_.r_cut;
    Real epsilon = params_.epsilon;
    Real sigma = params_.sigma;
    int power = params_.power;

    for (int a = 0; a < target.numNodes(); ++a) {
        for (int b = 0; b < source.numNodes(); ++b) {
            Vec2d diff = pos_t[a] - pos_s[b];
            diff.x = minImageDx(diff.x);
            Real r = diff.norm();
            if (r < 1e-15 || r >= r_cut) continue;

            Real ratio = sigma / r;
            Real ratio_pow = std::pow(ratio, power);
            Real F_mag = epsilon * ratio_pow;

            Vec2d r_hat = diff / r;
            f_t[a] += r_hat * F_mag; // only write to target
        }
    }
}

void RepulsionForce::computeAll(CapsuleSystem& system, Real y_bottom, Real y_top) {
    int n = system.numCapsules();
    if (n == 0) return;

    // Use cell-list acceleration when there are enough capsules to benefit.
    // The cell list reduces O(N_nodes^2) pair checks to O(N_nodes * avg_neighbors).
    if (system.totalNodes() > 100) {
        computeAllCellList(system, y_bottom, y_top);
        return;
    }

#ifdef _OPENMP
    // Parallel: one-sided force — each thread owns its target capsule.
    // Doubles FLOPs vs Newton-3rd-law but eliminates write conflicts.
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            computeOneSidedRepulsion(system[i], system[j]);
        }
        computeWallRepulsion(system[i], y_bottom, y_top);
    }
#else
    // Serial: use Newton's 3rd law (compute each pair once)
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            computeInterCapsule(system[i], system[j]);
        }
    }
    for (int i = 0; i < n; ++i) {
        computeWallRepulsion(system[i], y_bottom, y_top);
    }
#endif
}

void RepulsionForce::computeAllCellList(CapsuleSystem& system,
                                         Real y_bottom, Real y_top) {
    const int ncaps = system.numCapsules();
    const Real r_cut = params_.r_cut;
    const Real epsilon = params_.epsilon;
    const Real sigma = params_.sigma;
    const int power = params_.power;

    // Build the cell list from current node positions.
    // Cell size = r_cut so any interacting pair shares the same or adjacent cells.
    cell_list_.build(system, r_cut, Lx_, y_bottom, y_top);

    // One-sided, cell-list accelerated repulsion.
    // Outer loop: iterate over each capsule (and each of its nodes).
    // For each node, look up its cell and the 3x3 neighbor stencil,
    // then only check nodes from OTHER capsules found in those cells.
    //
    // Thread safety: each thread owns capsule i and only writes to
    // system[i].forces(), so no write conflicts occur.

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int i = 0; i < ncaps; ++i) {
        auto& cap_i = system[i];
        auto& pos_i = cap_i.positions();
        auto& f_i = cap_i.forces();
        const int nnodes_i = cap_i.numNodes();

        for (int a = 0; a < nnodes_i; ++a) {
            int ci = cell_list_.cellIndex(pos_i[a]);
            int cx = ci % cell_list_.ncx();
            int cy = ci / cell_list_.ncx();

            int neighbors[9];
            int nneigh = cell_list_.neighborCells(cx, cy, neighbors);

            for (int nc = 0; nc < nneigh; ++nc) {
                const auto& cell_nodes = cell_list_.cell(
                    neighbors[nc] % cell_list_.ncx(),
                    neighbors[nc] / cell_list_.ncx());

                for (const auto& ref : cell_nodes) {
                    if (ref.capsule == i) continue; // skip self-capsule

                    const Vec2d& pos_b = system[ref.capsule].nodePosition(ref.node);
                    Vec2d diff = pos_i[a] - pos_b;
                    diff.x = minImageDx(diff.x);
                    Real r = diff.norm();

                    if (r < 1e-15 || r >= r_cut) continue;

                    Real ratio = sigma / r;
                    Real ratio_pow = std::pow(ratio, power);
                    Real F_mag = epsilon * ratio_pow;

                    Vec2d r_hat = diff / r;
                    f_i[a] += r_hat * F_mag; // one-sided: only write to target
                }
            }
        }

        // Wall repulsion (always needed, independent of cell list)
        computeWallRepulsion(cap_i, y_bottom, y_top);
    }
}

void RepulsionForce::computeObstacleRepulsion(CapsuleSystem& system, const Obstacle& obs) {
    Real r_cut = params_.r_cut;
    Real epsilon = params_.epsilon;
    Real sigma = params_.sigma;
    int power = params_.power;
    int ncaps = system.numCapsules();

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic)
#endif
    for (int c = 0; c < ncaps; ++c) {
        auto& pos = system[c].positions();
        auto& f = system[c].forces();

        for (int i = 0; i < system[c].numNodes(); ++i) {
            Real sd = obs.signedDistance(pos[i].x, pos[i].y);

            // sd > 0: outside obstacle (normal repulsion)
            // sd < 0: INSIDE obstacle (emergency push-out)
            if (sd >= r_cut) continue;  // too far away

            Real r = std::abs(sd);
            if (r < 1e-15) r = 1e-15;  // avoid division by zero

            Vec2d normal = obs.normalAt(pos[i].x, pos[i].y);
            Real ratio = sigma / r;
            Real ratio_pow = std::pow(ratio, power);
            Real F_mag = epsilon * ratio_pow;

            // If inside obstacle, apply stronger emergency push-out
            if (sd < 0.0) {
                F_mag *= 10.0;
            }

            f[i] += normal * F_mag;
        }
    }
}

} // namespace softflow
