#pragma once

#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "lsmps/lsmps_matrices.hpp"
#include "neighbor/neighbor_search.hpp"
#include "provisional/provisional.hpp"

#include <vector>

namespace lsmps {

struct PressurePoissonSolveInfo {
    bool converged = false;
    int iterations = 0;
    double final_residual_norm = 0.0;
    int converged_reason = 0;
};

struct PressurePoissonDiagnostics {
    std::vector<double> rhs;
    std::vector<double> divergence;
    std::vector<double> wall_neumann_source;
    std::vector<int> pressure_dof;
};

struct PressurePoissonResult {
    bool solved = false;
    std::vector<double> pressure;
    PressurePoissonSolveInfo solve_info;
    PressurePoissonDiagnostics diagnostics;
};

class PressurePoissonAssembler {
public:
    PressurePoissonResult solve(
        const ParticleSet& particles,
        const TypedNeighborList& neighbors,
        const LsmpsMatrixSet& matrices,
        const ProvisionalVelocityResult& provisional,
        const SimulationConfig& config) const;
};

}  // namespace lsmps
