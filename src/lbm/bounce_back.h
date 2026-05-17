#pragma once
#include "lattice_field.h"
#include <vector>

namespace softflow {

class BounceBack {
public:
    // Apply halfway bounce-back for all SOLID cells.
    // Must be called after streaming (on the post-stream distributions).

    // Original method: iterates all nodes checking type (fallback)
    void apply(LatticeField& field);

    // Optimized method: uses precomputed list of solid node indices
    // Avoids per-step type-checking over all nodes.
    void apply(LatticeField& field, const std::vector<int>& solid_nodes);
};

} // namespace softflow
