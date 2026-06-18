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
    if (free_surface.screen_radius_factor <= 0.0) {
        errors.push_back("free_surface.screen_radius_factor must be positive");
    }
    if (free_surface.wall_patch_radius_factor <= 0.0) {
        errors.push_back("free_surface.wall_patch_radius_factor must be positive");
    }
    if (free_surface.particle_radius_factor <= 0.0) {
        errors.push_back("free_surface.particle_radius_factor must be positive");
    }
    if (free_surface.open_threshold < 0.0 || free_surface.open_threshold > 1.0) {
        errors.push_back("free_surface.open_threshold must be in [0, 1]");
    }
    if (free_surface.cone_angle_degrees <= 0.0 || free_surface.cone_angle_degrees >= 180.0) {
        errors.push_back("free_surface.cone_angle_degrees must be in (0, 180)");
    }
    if (free_surface.cone_threshold < 0.0 || free_surface.cone_threshold > 1.0) {
        errors.push_back("free_surface.cone_threshold must be in [0, 1]");
    }
    if (free_surface.min_cone_accessible_ratio < 0.0 || free_surface.min_cone_accessible_ratio > 1.0) {
        errors.push_back("free_surface.min_cone_accessible_ratio must be in [0, 1]");
    }
    if (free_surface.cubed_sphere_q == 0) {
        errors.push_back("free_surface.cubed_sphere_q must be positive");
    }
    if (free_surface.splash_open_threshold < 0.0 || free_surface.splash_open_threshold > 1.0) {
        errors.push_back("free_surface.splash_open_threshold must be in [0, 1]");
    }
    if (free_surface.near_surface_distance_factor <= 0.0) {
        errors.push_back("free_surface.near_surface_distance_factor must be positive");
    }

    if (lsmps.min_neighbors == 0) {
        errors.push_back("lsmps.min_neighbors must be positive");
    }
    if (lsmps.eigenvalue_tolerance <= 0.0) {
        errors.push_back("lsmps.eigenvalue_tolerance must be positive");
    }
    if (lsmps.condition_number_warning <= 0.0) {
        errors.push_back("lsmps.condition_number_warning must be positive");
    }
    if (lsmps.condition_number_failure <= 0.0) {
        errors.push_back("lsmps.condition_number_failure must be positive");
    }
    if (lsmps.condition_number_warning > lsmps.condition_number_failure) {
        errors.push_back("lsmps.condition_number_warning must not exceed condition_number_failure");
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
