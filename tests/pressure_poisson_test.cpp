#include "core/particle_set.hpp"
#include "core/particle_types.hpp"
#include "core/simulation_config.hpp"
#include "core/vector3.hpp"
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

constexpr int grid_size = 8;
constexpr double spacing = 1.0;
constexpr double support_radius = 3.1 * spacing;
constexpr double density = 1000.0;
constexpr double gravity = 9.81;

void addFluidBlock(lsmps::ParticleSet& particles, int nx, int ny, int nz) {
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const auto state = (y == ny - 1) ? lsmps::FluidParticleState::FreeSurface : lsmps::FluidParticleState::Internal;
                particles.addFluidParticle({spacing * x, spacing * y, spacing * z}, {}, state);
            }
        }
    }
}

void addBoxWalls(lsmps::ParticleSet& particles, int nx, int ny, int nz) {
    for (int z = -1; z <= nz; ++z) {
        for (int y = -1; y <= ny; ++y) {
            for (int x = -1; x <= nx; ++x) {
                const bool in_fluid = x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz;
                if (in_fluid || y == ny) {
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

lsmps::SimulationConfig testConfig() {
    lsmps::SimulationConfig config;
    config.time.dt = 0.01;
    config.geometry.particle_spacing = spacing;
    config.geometry.support_radius = support_radius;
    config.physical.density = density;
    config.physical.viscosity = 0.0;
    config.physical.gravity = {0.0, -gravity, 0.0};
    config.lsmps.min_neighbors = 9;
    config.lsmps.eigenvalue_tolerance = 1.0e-12;
    config.lsmps.condition_number_warning = 1.0e10;
    config.lsmps.condition_number_failure = 1.0e14;
    config.linear_solver.max_iterations = 1000;
    config.linear_solver.tolerance = 1.0e-10;
    return config;
}

std::vector<double> pressureError(const lsmps::ParticleSet& particles, const std::vector<double>& pressure) {
    std::vector<double> error(particles.size(), 0.0);
    const double height = static_cast<double>(grid_size - 1) * spacing;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            const double exact = density * gravity * (height - particles.positions()[i].y);
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

std::vector<double> intFieldAsDouble(const std::vector<int>& values) {
    std::vector<double> result(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        result[i] = static_cast<double>(values[i]);
    }
    return result;
}

}  // namespace

int main() {
    lsmps::ParticleSet particles;
    particles.reserve(1200);
    addFluidBlock(particles, grid_size, grid_size, grid_size);
    addBoxWalls(particles, grid_size, grid_size, grid_size);

    const lsmps::SimulationConfig config = testConfig();
    lsmps::NeighborSearch search(config.geometry.support_radius, {-spacing, -spacing, -spacing});
    const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
    search.updateNeighborCounts(particles, neighbors);
    const lsmps::LsmpsMatrixSet matrices =
        lsmps::buildLsmpsMatrices(particles, neighbors, config.geometry.support_radius, config.lsmps);

    const lsmps::ProvisionalVelocityCalculator provisional_calculator;
    const lsmps::ProvisionalVelocityResult provisional =
        provisional_calculator.compute(particles, neighbors, matrices, config);

    const lsmps::PressurePoissonAssembler assembler;
    const lsmps::PressurePoissonResult result =
        assembler.solve(particles, neighbors, matrices, provisional, config);
    assert(result.solved);
    assert(result.solve_info.converged);

    const std::vector<double> error = pressureError(particles, result.pressure);
    const double max_error = maxFluidError(particles, error);
    if (max_error >= 1.0e-5) {
        std::cerr << "PPE converged=" << result.solve_info.converged
                  << " reason=" << result.solve_info.converged_reason
                  << " iterations=" << result.solve_info.iterations
                  << " residual=" << result.solve_info.final_residual_norm
                  << " max_pressure_error=" << max_error << '\n';
    }
    assert(max_error < 1.0e-5);

    lsmps::ParticleSet output_particles = particles;
    output_particles.pressures() = result.pressure;
    const lsmps::VtkWriter writer;
    writer.writeParticles(
        "output/pressure_poisson_hydrostatic_8x8x8.vtk",
        output_particles,
        {
            {"pressure_error", error},
            {"ppe_rhs", result.diagnostics.rhs},
            {"ppe_divergence", result.diagnostics.divergence},
            {"ppe_wall_neumann_source", result.diagnostics.wall_neumann_source},
            {"pressure_dof", intFieldAsDouble(result.diagnostics.pressure_dof)},
        });

    return 0;
}
