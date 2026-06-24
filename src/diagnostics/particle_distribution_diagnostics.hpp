#pragma once

#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "neighbor/neighbor_search.hpp"

#include <cstddef>
#include <vector>

namespace lsmps {

struct ParticleDistributionDiagnostics {
    std::vector<double> nearest_fluid_distance;
    std::vector<double> max_fluid_neighbor_distance;
    std::vector<double> mean_fluid_neighbor_distance;
    std::vector<double> number_density;
    std::vector<double> number_density_ratio;
    std::vector<double> x_positive_neighbor_count;
    std::vector<double> x_negative_neighbor_count;
    std::vector<double> y_positive_neighbor_count;
    std::vector<double> y_negative_neighbor_count;
    std::vector<double> z_positive_neighbor_count;
    std::vector<double> z_negative_neighbor_count;
    std::vector<double> directional_coverage_score;
    std::vector<double> geometry_min_eigenvalue;
    std::vector<double> geometry_max_eigenvalue;
    std::vector<double> geometry_condition_number;
};

ParticleDistributionDiagnostics computeParticleDistributionDiagnostics(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    const GeometryConfig& geometry,
    LsmpsKernelType kernel_type);

}  // namespace lsmps
