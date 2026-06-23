#pragma once

#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "core/vector3.hpp"
#include "lsmps/lsmps_matrices.hpp"
#include "neighbor/neighbor_search.hpp"
#include "pressure_poisson/pressure_poisson.hpp"
#include "provisional/provisional.hpp"

#include <vector>

namespace lsmps {

enum class CorrectionParticleStatus {
    Updated = 0,
    MatrixUnavailable = 1,
    NonFluid = 2,
};

struct CorrectionResult {
    bool applied = false;
    std::vector<Vector3> pressure_gradient;
    std::vector<Vector3> velocity_correction;
    std::vector<Vector3> next_velocity;
    std::vector<Vector3> displacement;
    std::vector<Vector3> next_position;
    std::vector<CorrectionParticleStatus> status;
};

class PressureCorrectionApplier {
public:
    CorrectionResult apply(
        const ParticleSet& particles,
        const TypedNeighborList& neighbors,
        const LsmpsMatrixSet& matrices,
        const ProvisionalVelocityResult& provisional,
        const PressurePoissonResult& pressure,
        const SimulationConfig& config) const;
};

}  // namespace lsmps
