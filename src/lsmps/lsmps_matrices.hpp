#pragma once

#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "lsmps/lsmps_basis.hpp"
#include "neighbor/neighbor_search.hpp"

#include <cstddef>
#include <vector>

namespace lsmps {

enum class LsmpsMatrixStatus {
    Valid = 0,
    NotEnoughNeighbors = 1,
    RankDeficient = 2,
    IllConditioned = 3,
    InversionFailed = 4,
};

struct LsmpsInverseMatrix {
    LsmpsMomentMatrix inverse_moment = LsmpsMomentMatrix::Zero();
    LsmpsMatrixStatus status = LsmpsMatrixStatus::NotEnoughNeighbors;
    int rank = 0;
    double condition_number = 0.0;
    double min_eigenvalue = 0.0;
    double max_eigenvalue = 0.0;
    std::size_t total_neighbor_count = 0;
    std::size_t fluid_neighbor_count = 0;
    std::size_t wall_neighbor_count = 0;
};

struct LsmpsParticleMatrices {
    LsmpsInverseMatrix regular;
    LsmpsInverseMatrix pressure_neumann;
};

struct LsmpsMatrixSet {
    std::vector<LsmpsParticleMatrices> particles;
};

LsmpsMatrixSet buildLsmpsMatrices(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    double support_radius,
    const LsmpsConfig& config);

}  // namespace lsmps
