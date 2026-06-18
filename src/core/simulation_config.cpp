#include "core/simulation_config.hpp"

#include <string>
#include <vector>

namespace lsmps {
namespace {

bool isStrictlyLess(const Vector3& lhs, const Vector3& rhs) {
    return lhs.x < rhs.x && lhs.y < rhs.y && lhs.z < rhs.z;
}

}  // namespace

std::vector<std::string> SimulationConfig::validate() const {
    std::vector<std::string> errors;

    if (time.dt <= 0.0) {
        errors.push_back("time.dt must be positive");
    }
    if (time.end_time <= 0.0) {
        errors.push_back("time.end_time must be positive");
    }
    if (time.output_interval <= 0.0) {
        errors.push_back("time.output_interval must be positive");
    }
    if (time.dt > time.end_time) {
        errors.push_back("time.dt must not be larger than time.end_time");
    }

    if (geometry.particle_spacing <= 0.0) {
        errors.push_back("geometry.particle_spacing must be positive");
    }
    if (geometry.support_radius <= 0.0) {
        errors.push_back("geometry.support_radius must be positive");
    }
    if (geometry.support_radius < geometry.particle_spacing) {
        errors.push_back("geometry.support_radius must not be smaller than geometry.particle_spacing");
    }
    if (!isStrictlyLess(geometry.domain_min, geometry.domain_max)) {
        errors.push_back("geometry.domain_min must be smaller than geometry.domain_max in every component");
    }
    if (!isStrictlyLess(geometry.fluid_min, geometry.fluid_max)) {
        errors.push_back("geometry.fluid_min must be smaller than geometry.fluid_max in every component");
    }

    if (physical.density <= 0.0) {
        errors.push_back("physical.density must be positive");
    }
    if (physical.viscosity < 0.0) {
        errors.push_back("physical.viscosity must be non-negative");
    }

    if (free_surface.neighbor_count_ratio <= 0.0 || free_surface.neighbor_count_ratio > 1.0) {
        errors.push_back("free_surface.neighbor_count_ratio must be in (0, 1]");
    }
    if (free_surface.number_density_ratio <= 0.0 || free_surface.number_density_ratio > 1.0) {
        errors.push_back("free_surface.number_density_ratio must be in (0, 1]");
    }

    if (linear_solver.max_iterations == 0) {
        errors.push_back("linear_solver.max_iterations must be positive");
    }
    if (linear_solver.tolerance <= 0.0) {
        errors.push_back("linear_solver.tolerance must be positive");
    }

    return errors;
}

}  // namespace lsmps
