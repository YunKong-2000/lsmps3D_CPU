#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "core/vector3.hpp"
#include "io/vtk_writer.hpp"
#include "lsmps/lsmps_matrices.hpp"
#include "neighbor/neighbor_search.hpp"
#include "provisional/provisional.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace {

constexpr int grid_size = 12;
constexpr double spacing = 1.0;
constexpr double support_radius = 3.1 * spacing;

void addFluidBlock(lsmps::ParticleSet& particles, int nx, int ny, int nz) {
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                particles.addFluidParticle({spacing * x, spacing * y, spacing * z});
            }
        }
    }
}

void addBoxWalls(lsmps::ParticleSet& particles, int nx, int ny, int nz) {
    for (int z = -1; z <= nz; ++z) {
        for (int y = -1; y <= ny; ++y) {
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
                if (y == ny) {
                    normal += {0.0, -1.0, 0.0};
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
                    particles.addWallParticle({spacing * x, spacing * y, spacing * z}, {}, 0.0, 0.0, normal);
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
    config.physical.viscosity = 0.25;
    config.physical.gravity = {0.0, -9.81, 0.0};
    config.lsmps.min_neighbors = 9;
    config.lsmps.eigenvalue_tolerance = 1.0e-12;
    config.lsmps.condition_number_warning = 1.0e10;
    config.lsmps.condition_number_failure = 1.0e14;
    return config;
}

std::vector<double> statusField(const lsmps::ProvisionalVelocityResult& result) {
    std::vector<double> values(result.status.size(), 0.0);
    for (std::size_t i = 0; i < result.status.size(); ++i) {
        values[i] = static_cast<double>(static_cast<int>(result.status[i]));
    }
    return values;
}

std::vector<double> vectorMagnitude(const std::vector<lsmps::Vector3>& values) {
    std::vector<double> result(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        result[i] = lsmps::norm(values[i]);
    }
    return result;
}

std::vector<double> viscousError(
    const lsmps::ParticleSet& particles,
    const lsmps::ProvisionalVelocityResult& result,
    const lsmps::Vector3& expected) {
    std::vector<double> values(particles.size(), 0.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            values[i] = lsmps::norm(result.viscous_acceleration[i] - expected);
        }
    }
    return values;
}

double maxFluidValue(const lsmps::ParticleSet& particles, const std::vector<double>& values) {
    double maximum = 0.0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            maximum = std::max(maximum, values[i]);
        }
    }
    return maximum;
}

}  // namespace

int main() {
    lsmps::ParticleSet particles;
    particles.reserve(4096);
    addFluidBlock(particles, grid_size, grid_size, grid_size);
    addBoxWalls(particles, grid_size, grid_size, grid_size);

    for (std::size_t i = 0; i < particles.size(); ++i) {
        const lsmps::Vector3& x = particles.positions()[i];
        particles.velocities()[i] = {
            x.y * x.y - x.z * x.z,
            2.0 * x.z * x.z - x.x * x.x,
            3.0 * x.x * x.x - x.y * x.y,
        };
    }

    const lsmps::SimulationConfig config = testConfig();
    lsmps::NeighborSearch search(config.geometry.support_radius, {-spacing, -spacing, -spacing});
    const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
    search.updateNeighborCounts(particles, neighbors);
    const lsmps::LsmpsMatrixSet matrices =
        lsmps::buildLsmpsMatrices(particles, neighbors, config.geometry.support_radius, config.lsmps);

    const lsmps::ProvisionalVelocityCalculator calculator;
    const lsmps::ProvisionalVelocityResult result =
        calculator.compute(particles, neighbors, matrices, config);
    assert(result.computed);

    const lsmps::Vector3 expected_viscous = {
        0.0,
        2.0 * config.physical.viscosity,
        4.0 * config.physical.viscosity,
    };
    const std::vector<double> viscosity_error = viscousError(particles, result, expected_viscous);
    assert(maxFluidValue(particles, viscosity_error) < 1.0e-8);

    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!particles.isFluid(i)) {
            continue;
        }
        assert(result.status[i] == lsmps::ProvisionalParticleStatus::Updated);
        const lsmps::Vector3 expected_delta =
            config.time.dt * (expected_viscous + config.physical.gravity);
        assert(lsmps::norm(result.velocity_delta[i] - expected_delta) < 1.0e-8);
        assert(lsmps::norm(result.provisional_velocity[i] - particles.velocities()[i] - expected_delta) < 1.0e-8);
    }

    const lsmps::VtkWriter writer;
    writer.writeParticles(
        "output/provisional_velocity_12x12x12.vtk",
        particles,
        {
            {"provisional_status", statusField(result)},
            {"provisional_velocity_magnitude", vectorMagnitude(result.provisional_velocity)},
            {"viscous_acceleration_magnitude", vectorMagnitude(result.viscous_acceleration)},
            {"body_acceleration_magnitude", vectorMagnitude(result.body_acceleration)},
            {"velocity_delta_magnitude", vectorMagnitude(result.velocity_delta)},
            {"viscous_acceleration_error", viscosity_error},
        },
        {
            {"provisional_velocity", result.provisional_velocity},
            {"viscous_acceleration", result.viscous_acceleration},
            {"body_acceleration", result.body_acceleration},
            {"velocity_delta", result.velocity_delta},
        });

    return 0;
}
