#pragma once

#include <cstddef>

namespace lsmps {

struct SimulationState {
    std::size_t current_step = 0;
    double current_time = 0.0;
    double max_velocity = 0.0;
    double cfl_number = 0.0;
    double pressure_residual = 0.0;
    std::size_t min_neighbor_count = 0;
    std::size_t max_neighbor_count = 0;
    double mean_neighbor_count = 0.0;
};

}  // namespace lsmps
