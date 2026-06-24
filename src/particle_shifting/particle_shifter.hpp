#pragma once

#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "core/vector3.hpp"
#include "neighbor/neighbor_search.hpp"

#include <vector>

namespace lsmps {

struct ParticleShiftResult {
    bool applied = false;
    std::vector<Vector3> displacement;
    std::vector<double> magnitude;
    std::vector<double> limited;
    std::vector<double> repulsion_active;
};

class ParticleShifter {
public:
    ParticleShiftResult compute(
        const ParticleSet& particles,
        const TypedNeighborList& neighbors,
        const SimulationConfig& config) const;

    void apply(ParticleSet& particles, const ParticleShiftResult& result) const;
};

}  // namespace lsmps
