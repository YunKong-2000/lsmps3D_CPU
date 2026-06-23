#include "core/particle_set.hpp"
#include "core/particle_types.hpp"
#include "core/simulation_config.hpp"
#include "core/vector3.hpp"
#include "time_integration/time_stepper.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace {

constexpr int grid_size = 6;
constexpr double spacing = 1.0;
constexpr double density = 1000.0;
constexpr double gravity = 9.81;

lsmps::SimulationConfig testConfig() {
    lsmps::SimulationConfig config;
    config.time.dt = 0.005;
    config.time.end_time = 0.05;
    config.time.output_interval = 0.01;
    config.time_step.start_time = 0.0;
    config.time_step.end_time = config.time.end_time;
    config.time_step.initial_dt = config.time.dt;
    config.time_step.min_dt = 0.001;
    config.time_step.max_dt = config.time.dt;
    config.time_step.cfl_number = 0.2;
    config.time_step.growth_factor = 1.05;
    config.time_step.output_interval = config.time.output_interval;
    config.file.output_directory = "output/time_stepper_hydrostatic_long";
    config.file.output_prefix = "hydrostatic";
    config.geometry.particle_spacing = spacing;
    config.geometry.support_radius = 3.1 * spacing;
    config.geometry.domain_min = {-spacing, -spacing, -spacing};
    config.geometry.domain_max = {grid_size * spacing, grid_size * spacing, grid_size * spacing};
    config.physical.density = density;
    config.physical.viscosity = 0.0;
    config.physical.gravity = {0.0, -gravity, 0.0};
    config.free_surface.screen_radius_factor = 3.1;
    config.free_surface.wall_patch_radius_factor = 0.85;
    config.free_surface.particle_radius_factor = 0.5;
    config.free_surface.open_threshold = 0.18;
    config.free_surface.cone_angle_degrees = 45.0;
    config.free_surface.cone_threshold = 0.62;
    config.free_surface.min_cone_accessible_ratio = 0.40;
    config.free_surface.cubed_sphere_q = 8;
    config.free_surface.splash_max_fluid_neighbors = 4;
    config.free_surface.splash_open_threshold = 0.75;
    config.free_surface.near_surface_distance_factor = 1.5;
    config.lsmps.min_neighbors = 9;
    config.lsmps.eigenvalue_tolerance = 1.0e-12;
    config.lsmps.condition_number_warning = 1.0e10;
    config.lsmps.condition_number_failure = 1.0e14;
    config.linear_solver.max_iterations = 1000;
    config.linear_solver.tolerance = 1.0e-10;
    return config;
}

void addFluidBlock(lsmps::ParticleSet& particles) {
    for (int z = 0; z < grid_size; ++z) {
        for (int y = 0; y < grid_size; ++y) {
            for (int x = 0; x < grid_size; ++x) {
                particles.addFluidParticle(
                    {spacing * x, spacing * y, spacing * z},
                    {},
                    lsmps::FluidParticleState::Internal,
                    0.0,
                    density);
            }
        }
    }
}

void addBoxWalls(lsmps::ParticleSet& particles) {
    for (int z = -1; z <= grid_size; ++z) {
        for (int y = -1; y <= grid_size; ++y) {
            for (int x = -1; x <= grid_size; ++x) {
                const bool in_fluid =
                    x >= 0 && x < grid_size && y >= 0 && y < grid_size && z >= 0 && z < grid_size;
                if (in_fluid || y == grid_size) {
                    continue;
                }

                lsmps::Vector3 normal;
                bool is_wall = false;
                if (x == -1) {
                    normal += {1.0, 0.0, 0.0};
                    is_wall = true;
                }
                if (x == grid_size) {
                    normal += {-1.0, 0.0, 0.0};
                    is_wall = true;
                }
                if (y == -1) {
                    normal += {0.0, 1.0, 0.0};
                    is_wall = true;
                }
                if (z == -1) {
                    normal += {0.0, 0.0, 1.0};
                    is_wall = true;
                }
                if (z == grid_size) {
                    normal += {0.0, 0.0, -1.0};
                    is_wall = true;
                }

                if (is_wall) {
                    const double length = lsmps::norm(normal);
                    if (length > 0.0) {
                        normal /= length;
                    }
                    particles.addWallParticle({spacing * x, spacing * y, spacing * z}, {}, 0.0, density, normal);
                }
            }
        }
    }
}

double maxFluidSpeed(const lsmps::ParticleSet& particles) {
    double maximum = 0.0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            maximum = std::max(maximum, lsmps::norm(particles.velocities()[i]));
        }
    }
    return maximum;
}

double maxFluidPositionChange(const lsmps::ParticleSet& particles, const std::vector<lsmps::Vector3>& initial) {
    double maximum = 0.0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            maximum = std::max(maximum, lsmps::norm(particles.positions()[i] - initial[i]));
        }
    }
    return maximum;
}

}  // namespace

int main() {
    const lsmps::SimulationConfig config = testConfig();

    lsmps::ParticleSet particles;
    particles.reserve(700);
    addFluidBlock(particles);
    addBoxWalls(particles);
    const std::vector<lsmps::Vector3> initial_positions = particles.positions();

    lsmps::TimeStepper stepper(config);
    const std::vector<lsmps::TimeStepDiagnostics> history = stepper.run(particles);

    assert(history.size() == 10);
    for (const lsmps::TimeStepDiagnostics& diagnostics : history) {
        if (!diagnostics.pressure_converged || diagnostics.max_velocity >= 1.0e-6 ||
            diagnostics.cfl_number >= 1.0e-6 || diagnostics.free_surface_count == 0) {
            std::cerr << "step=" << diagnostics.step
                      << " pressure_converged=" << diagnostics.pressure_converged
                      << " iterations=" << diagnostics.pressure_iterations
                      << " residual=" << diagnostics.pressure_residual
                      << " max_velocity=" << diagnostics.max_velocity
                      << " cfl=" << diagnostics.cfl_number
                      << " free_surface_count=" << diagnostics.free_surface_count << '\n';
        }
        assert(diagnostics.pressure_converged);
        assert(diagnostics.max_velocity < 1.0e-6);
        assert(diagnostics.cfl_number < 1.0e-6);
        assert(diagnostics.free_surface_count > 0);
    }

    assert(maxFluidSpeed(particles) < 1.0e-6);
    assert(maxFluidPositionChange(particles, initial_positions) < 1.0e-8);

    return 0;
}
