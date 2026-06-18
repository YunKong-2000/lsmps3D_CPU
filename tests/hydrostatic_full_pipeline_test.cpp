#include "core/particle_set.hpp"
#include "core/particle_types.hpp"
#include "core/simulation_config.hpp"
#include "core/vector3.hpp"
#include "free_surface/free_surface_detector.hpp"
#include "io/vtk_writer.hpp"
#include "lsmps/lsmps_matrices.hpp"
#include "neighbor/neighbor_search.hpp"
#include "pressure_poisson/pressure_poisson.hpp"
#include "provisional/provisional.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

namespace {

constexpr double tank_size = 1.0;
constexpr double water_height = 0.5;
constexpr double spacing = 0.05;
constexpr int nx = static_cast<int>(tank_size / spacing) + 1;
constexpr int ny = static_cast<int>(water_height / spacing) + 1;
constexpr int nz = static_cast<int>(tank_size / spacing) + 1;
constexpr int tank_wall_ny = static_cast<int>(tank_size / spacing) + 1;
constexpr double density = 1000.0;
constexpr double gravity = 9.81;

lsmps::SimulationConfig hydrostaticConfig() {
    lsmps::SimulationConfig config;
    config.time.dt = 0.001;
    config.geometry.particle_spacing = spacing;
    config.geometry.support_radius = 3.1 * spacing;
    config.geometry.domain_min = {0.0, 0.0, 0.0};
    config.geometry.domain_max = {tank_size, tank_size, tank_size};
    config.geometry.fluid_min = {0.0, 0.0, 0.0};
    config.geometry.fluid_max = {tank_size, water_height, tank_size};
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
    config.linear_solver.max_iterations = 2000;
    config.linear_solver.tolerance = 1.0e-10;
    return config;
}

void addFluid(lsmps::ParticleSet& particles) {
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                particles.addFluidParticle({spacing * x, spacing * y, spacing * z}, {}, lsmps::FluidParticleState::Internal, 0.0, density);
            }
        }
    }
}

void addTankWalls(lsmps::ParticleSet& particles) {
    for (int z = -1; z <= nz; ++z) {
        for (int y = -1; y < tank_wall_ny; ++y) {
            for (int x = -1; x <= nx; ++x) {
                const bool in_fluid = x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz;
                if (in_fluid) {
                    continue;
                }

                lsmps::Vector3 normal;
                bool is_wall = false;
                if (x == -1) {
                    normal += {1.0, 0.0, 0.0};
                    is_wall = true;
                }
                if (x == nx) {
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
                if (z == nz) {
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

std::vector<double> pressureError(const lsmps::ParticleSet& particles, const std::vector<double>& pressure) {
    std::vector<double> error(particles.size(), 0.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            const double exact = density * gravity * (water_height - particles.positions()[i].y);
            error[i] = std::abs(pressure[i] - exact);
        }
    }
    return error;
}

double maxFluidError(const lsmps::ParticleSet& particles, const std::vector<double>& error) {
    double result = 0.0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            result = std::max(result, error[i]);
        }
    }
    return result;
}

std::size_t countState(const lsmps::ParticleSet& particles, lsmps::FluidParticleState state) {
    std::size_t count = 0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i) && particles.fluidStates()[i] == state) {
            ++count;
        }
    }
    return count;
}

std::vector<double> intFieldAsDouble(const std::vector<int>& values) {
    std::vector<double> result(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        result[i] = static_cast<double>(values[i]);
    }
    return result;
}

std::vector<double> stateAsDouble(const lsmps::ParticleSet& particles) {
    std::vector<double> result(particles.size(), 0.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        result[i] = static_cast<double>(static_cast<int>(particles.fluidStates()[i]));
    }
    return result;
}

}  // namespace

int main() {
    const lsmps::SimulationConfig config = hydrostaticConfig();

    lsmps::ParticleSet particles;
    particles.reserve(8000);
    addFluid(particles);
    addTankWalls(particles);

    lsmps::NeighborSearch search(config.geometry.support_radius, {-spacing, -spacing, -spacing});
    const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
    search.updateNeighborCounts(particles, neighbors);

    const lsmps::FreeSurfaceDetector detector(config.free_surface, config.geometry.particle_spacing);
    const lsmps::FreeSurfaceDiagnostics free_surface = detector.detect(particles, neighbors);

    assert(countState(particles, lsmps::FluidParticleState::FreeSurface) > 0);
    assert(countState(particles, lsmps::FluidParticleState::NearFreeSurface) > 0);

    const lsmps::LsmpsMatrixSet matrices =
        lsmps::buildLsmpsMatrices(particles, neighbors, config.geometry.support_radius, config.lsmps);

    const lsmps::ProvisionalVelocityCalculator provisional_calculator;
    const lsmps::ProvisionalVelocityResult provisional =
        provisional_calculator.compute(particles, neighbors, matrices, config);
    assert(provisional.computed);

    const lsmps::PressurePoissonAssembler pressure_poisson;
    const lsmps::PressurePoissonResult pressure =
        pressure_poisson.solve(particles, neighbors, matrices, provisional, config);
    assert(pressure.solved);
    assert(pressure.solve_info.converged);

    const std::vector<double> error = pressureError(particles, pressure.pressure);
    const double max_error = maxFluidError(particles, error);
    if (max_error > 1.0e-4) {
        std::cerr << "hydrostatic full pipeline max pressure error = " << max_error
                  << ", iterations = " << pressure.solve_info.iterations
                  << ", residual = " << pressure.solve_info.final_residual_norm << '\n';
    }
    assert(max_error < 1.0e-4);

    lsmps::ParticleSet output_particles = particles;
    output_particles.pressures() = pressure.pressure;
    const lsmps::VtkWriter writer;
    writer.writeParticles(
        "output/hydrostatic_box_1m_0p5m_full_pipeline.vtk",
        output_particles,
        {
            {"pressure_error", error},
            {"fluid_state_explicit", stateAsDouble(output_particles)},
            {"free_surface_open_ratio", free_surface.open_ratio},
            {"free_surface_cone_ratio", free_surface.cone_ratio},
            {"free_surface_reason_code", intFieldAsDouble(free_surface.reason_code)},
            {"ppe_rhs", pressure.diagnostics.rhs},
            {"ppe_divergence", pressure.diagnostics.divergence},
            {"ppe_wall_neumann_source", pressure.diagnostics.wall_neumann_source},
            {"pressure_dof", intFieldAsDouble(pressure.diagnostics.pressure_dof)},
        },
        {
            {"provisional_velocity", provisional.provisional_velocity},
            {"velocity_delta", provisional.velocity_delta},
        });

    return 0;
}
