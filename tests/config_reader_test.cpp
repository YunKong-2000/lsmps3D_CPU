#include "config/config_reader.hpp"
#include "core/simulation_config.hpp"

#include <cassert>
#include <stdexcept>
#include <string>

int main() {
    const lsmps::SimulationConfig defaults;
    assert(defaults.time.dt == 0.001);
    assert(defaults.geometry.particle_spacing == 0.02);
    assert(defaults.geometry.support_radius == 0.06);
    assert(defaults.physical.density == 1000.0);
    assert(defaults.physical.gravity.y == -9.81);
    assert(defaults.validate().empty());

    const lsmps::ConfigReader reader;
    const lsmps::SimulationConfig config = reader.readString(R"(
        # Time parameters
        [time]
        dt = 0.005
        end_time = 2.5
        output_interval = 0.1

        [geometry]
        particle_spacing = 0.025
        support_radius = 0.075
        domain_min = -1.0 0.0 -2.0
        domain_max = 1.0 2.0 2.0
        fluid_min = -0.5 0.0 -0.5
        fluid_max = 0.5 1.0 0.5

        [physical]
        density = 998.2
        viscosity = 0.0000011
        gravity = 0.0 -9.8 0.0

        [free_surface]
        neighbor_count_ratio = 0.8
        number_density_ratio = 0.75
        near_surface_layers = 2
        screen_radius_factor = 3.1
        wall_patch_radius_factor = 0.85
        particle_radius_factor = 0.5
        open_threshold = 0.18
        cone_angle_degrees = 45.0
        cone_threshold = 0.62
        min_cone_accessible_ratio = 0.40
        cubed_sphere_q = 8
        splash_max_fluid_neighbors = 3
        splash_open_threshold = 0.8
        near_surface_distance_factor = 1.25

        [lsmps]
        min_neighbors = 12
        eigenvalue_tolerance = 0.000000000001
        condition_number_warning = 1000000
        condition_number_failure = 1000000000
        kernel_type = linear

        [linear_solver]
        max_iterations = 500
        tolerance = 0.000001
    )");

    assert(config.time.dt == 0.005);
    assert(config.time.end_time == 2.5);
    assert(config.time.output_interval == 0.1);
    assert(config.geometry.particle_spacing == 0.025);
    assert(config.geometry.support_radius == 0.075);
    assert(config.geometry.domain_min.x == -1.0);
    assert(config.geometry.domain_max.z == 2.0);
    assert(config.geometry.fluid_max.y == 1.0);
    assert(config.physical.density == 998.2);
    assert(config.physical.viscosity == 0.0000011);
    assert(config.physical.gravity.y == -9.8);
    assert(config.free_surface.neighbor_count_ratio == 0.8);
    assert(config.free_surface.number_density_ratio == 0.75);
    assert(config.free_surface.near_surface_layers == 2);
    assert(config.free_surface.screen_radius_factor == 3.1);
    assert(config.free_surface.wall_patch_radius_factor == 0.85);
    assert(config.free_surface.particle_radius_factor == 0.5);
    assert(config.free_surface.open_threshold == 0.18);
    assert(config.free_surface.cone_angle_degrees == 45.0);
    assert(config.free_surface.cone_threshold == 0.62);
    assert(config.free_surface.min_cone_accessible_ratio == 0.40);
    assert(config.free_surface.cubed_sphere_q == 8);
    assert(config.free_surface.splash_max_fluid_neighbors == 3);
    assert(config.free_surface.splash_open_threshold == 0.8);
    assert(config.free_surface.near_surface_distance_factor == 1.25);
    assert(config.lsmps.min_neighbors == 12);
    assert(config.lsmps.eigenvalue_tolerance == 0.000000000001);
    assert(config.lsmps.condition_number_warning == 1000000.0);
    assert(config.lsmps.condition_number_failure == 1000000000.0);
    assert(config.lsmps.kernel_type == lsmps::LsmpsKernelType::Linear);
    assert(config.linear_solver.max_iterations == 500);
    assert(config.linear_solver.tolerance == 0.000001);

    bool rejected_invalid_value = false;
    try {
        reader.readString("[time]\ndt = -1.0\n");
    } catch (const std::runtime_error&) {
        rejected_invalid_value = true;
    }
    assert(rejected_invalid_value);

    bool rejected_unknown_key = false;
    try {
        reader.readString("[time]\nunknown = 1.0\n");
    } catch (const std::runtime_error&) {
        rejected_unknown_key = true;
    }
    assert(rejected_unknown_key);

    bool rejected_negative_integer = false;
    try {
        reader.readString("[free_surface]\nnear_surface_layers = -1\n");
    } catch (const std::runtime_error&) {
        rejected_negative_integer = true;
    }
    assert(rejected_negative_integer);

    bool rejected_missing_section = false;
    try {
        reader.readString("dt = 0.001\n");
    } catch (const std::runtime_error&) {
        rejected_missing_section = true;
    }
    assert(rejected_missing_section);

    bool rejected_prefixed_key = false;
    try {
        reader.readString("[time]\ntime.dt = 0.001\n");
    } catch (const std::runtime_error&) {
        rejected_prefixed_key = true;
    }
    assert(rejected_prefixed_key);

    return 0;
}
