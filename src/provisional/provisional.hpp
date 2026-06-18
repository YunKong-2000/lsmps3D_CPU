#pragma once

#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "core/vector3.hpp"
#include "lsmps/lsmps_matrices.hpp"
#include "neighbor/neighbor_search.hpp"

#include <vector>

namespace lsmps {

enum class ProvisionalParticleStatus {
    Updated = 0,
    MatrixUnavailable = 1,
    NonFluidOrWall = 2,
};

struct ProvisionalVelocityResult {
    bool computed = false;
    std::vector<Vector3> provisional_velocity;
    std::vector<Vector3> viscous_acceleration;
    std::vector<Vector3> body_acceleration;
    std::vector<Vector3> velocity_delta;
    std::vector<ProvisionalParticleStatus> status;
};

class ProvisionalVelocityCalculator {
public:
    ProvisionalVelocityResult compute(
        const ParticleSet& particles,
        const TypedNeighborList& neighbors,
        const LsmpsMatrixSet& matrices,
        const SimulationConfig& config,
        const std::vector<Vector3>& external_acceleration = {}) const;
};

}  // namespace lsmps
