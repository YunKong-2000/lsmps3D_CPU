#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "core/vector3.hpp"
#include "io/vtk_writer.hpp"
#include "lsmps/lsmps_basis.hpp"
#include "lsmps/lsmps_matrices.hpp"
#include "lsmps/weight_function.hpp"
#include "neighbor/neighbor_search.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr int grid_size = 20;
constexpr double spacing = 1.0;
constexpr double support_radius = 3.1 * spacing;
constexpr double density = 1000.0;
constexpr double gravity_magnitude = 9.81;

struct ScalarOperatorResult {
    std::vector<lsmps::Vector3> gradient;
    std::vector<double> laplacian;
    std::vector<double> gradient_error;
    std::vector<double> laplacian_error;
};

struct VectorOperatorResult {
    std::vector<double> grad_xx;
    std::vector<double> grad_xy;
    std::vector<double> grad_xz;
    std::vector<double> grad_yx;
    std::vector<double> grad_yy;
    std::vector<double> grad_yz;
    std::vector<double> grad_zx;
    std::vector<double> grad_zy;
    std::vector<double> grad_zz;
    std::vector<double> divergence;
    std::vector<double> laplacian_x;
    std::vector<double> laplacian_y;
    std::vector<double> laplacian_z;
    std::vector<double> divergence_error;
    std::vector<double> laplacian_error;
};

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

lsmps::LsmpsConfig validationConfig() {
    lsmps::LsmpsConfig config;
    config.min_neighbors = 9;
    config.eigenvalue_tolerance = 1.0e-12;
    config.condition_number_warning = 1.0e10;
    config.condition_number_failure = 1.0e14;
    return config;
}

std::vector<double> matrixStatusField(
    const lsmps::ParticleSet& particles,
    const lsmps::LsmpsMatrixSet& matrices,
    bool pressure) {
    std::vector<double> values(particles.size(), -1.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const auto status = pressure ? matrices.particles[i].pressure_neumann.status : matrices.particles[i].regular.status;
        values[i] = static_cast<double>(static_cast<int>(status));
    }
    return values;
}

std::vector<double> matrixConditionField(
    const lsmps::ParticleSet& particles,
    const lsmps::LsmpsMatrixSet& matrices,
    bool pressure) {
    std::vector<double> values(particles.size(), 0.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        values[i] = pressure ? matrices.particles[i].pressure_neumann.condition_number : matrices.particles[i].regular.condition_number;
    }
    return values;
}

std::vector<double> matrixRankField(
    const lsmps::ParticleSet& particles,
    const lsmps::LsmpsMatrixSet& matrices,
    bool pressure) {
    std::vector<double> values(particles.size(), 0.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        values[i] = static_cast<double>(pressure ? matrices.particles[i].pressure_neumann.rank : matrices.particles[i].regular.rank);
    }
    return values;
}

std::vector<double> wallCountField(
    const lsmps::ParticleSet& particles,
    const lsmps::LsmpsMatrixSet& matrices) {
    std::vector<double> values(particles.size(), 0.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        values[i] = static_cast<double>(matrices.particles[i].regular.wall_neighbor_count);
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

lsmps::LsmpsBasisVector scalarRhs(
    const lsmps::ParticleSet& particles,
    const lsmps::TypedNeighborList& neighbors,
    std::size_t i,
    const std::vector<double>& values,
    const lsmps::LsmpsConfig& config) {
    lsmps::LsmpsBasisVector rhs = lsmps::LsmpsBasisVector::Zero();
    const lsmps::Vector3 position_i = particles.positions()[i];

    for (const std::size_t j : neighbors.fluid[i]) {
        const lsmps::Vector3 offset = particles.positions()[j] - position_i;
        const double weight = lsmps::evaluateWeight(lsmps::norm(offset), support_radius, config.kernel_type);
        rhs.noalias() += weight * lsmps::evaluateTypeABasis(offset, support_radius) * (values[j] - values[i]);
    }
    for (const std::size_t j : neighbors.wall[i]) {
        const lsmps::Vector3 offset = particles.positions()[j] - position_i;
        const double weight = lsmps::evaluateWeight(lsmps::norm(offset), support_radius, config.kernel_type);
        rhs.noalias() += weight * lsmps::evaluateTypeABasis(offset, support_radius) * (values[j] - values[i]);
    }
    return rhs;
}

lsmps::LsmpsBasisVector pressureRhs(
    const lsmps::ParticleSet& particles,
    const lsmps::TypedNeighborList& neighbors,
    std::size_t i,
    const std::vector<double>& pressure,
    const lsmps::LsmpsConfig& config) {
    lsmps::LsmpsBasisVector rhs = lsmps::LsmpsBasisVector::Zero();
    const lsmps::Vector3 position_i = particles.positions()[i];

    for (const std::size_t j : neighbors.fluid[i]) {
        const lsmps::Vector3 offset = particles.positions()[j] - position_i;
        const double weight = lsmps::evaluateWeight(lsmps::norm(offset), support_radius, config.kernel_type);
        rhs.noalias() += weight * lsmps::evaluateTypeABasis(offset, support_radius) * (pressure[j] - pressure[i]);
    }
    for (const std::size_t j : neighbors.wall[i]) {
        const lsmps::Vector3 offset = particles.positions()[j] - position_i;
        const double weight = lsmps::evaluateWeight(lsmps::norm(offset), support_radius, config.kernel_type);
        const double neumann = lsmps::dot({0.0, -density * gravity_magnitude, 0.0}, particles.wallNormals()[j]);
        rhs.noalias() += weight *
                          lsmps::evaluateTypeANeumannBasis(offset, particles.wallNormals()[j], support_radius) *
                          support_radius * neumann;
    }
    return rhs;
}

lsmps::Vector3 gradientFromD(const lsmps::LsmpsBasisVector& d) {
    return {d[0] / support_radius, d[1] / support_radius, d[2] / support_radius};
}

double laplacianFromD(const lsmps::LsmpsBasisVector& d) {
    return (d[3] + d[4] + d[5]) / (support_radius * support_radius);
}

ScalarOperatorResult computePressureOperators(
    const lsmps::ParticleSet& particles,
    const lsmps::TypedNeighborList& neighbors,
    const lsmps::LsmpsMatrixSet& matrices,
    const std::vector<double>& pressure,
    const lsmps::LsmpsConfig& config) {
    ScalarOperatorResult result;
    result.gradient.assign(particles.size(), {});
    result.laplacian.assign(particles.size(), 0.0);
    result.gradient_error.assign(particles.size(), 0.0);
    result.laplacian_error.assign(particles.size(), 0.0);

    const lsmps::Vector3 analytical_gradient = {0.0, -density * gravity_magnitude, 0.0};
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!particles.isFluid(i) ||
            matrices.particles[i].pressure_neumann.status != lsmps::LsmpsMatrixStatus::Valid) {
            continue;
        }

        const lsmps::LsmpsBasisVector d =
            matrices.particles[i].pressure_neumann.inverse_moment * pressureRhs(particles, neighbors, i, pressure, config);
        result.gradient[i] = gradientFromD(d);
        result.laplacian[i] = laplacianFromD(d);
        result.gradient_error[i] = lsmps::norm(result.gradient[i] - analytical_gradient);
        result.laplacian_error[i] = std::abs(result.laplacian[i]);
    }
    return result;
}

VectorOperatorResult computeVelocityOperators(
    const lsmps::ParticleSet& particles,
    const lsmps::TypedNeighborList& neighbors,
    const lsmps::LsmpsMatrixSet& matrices,
    const std::vector<double>& ux,
    const std::vector<double>& uy,
    const std::vector<double>& uz,
    const lsmps::LsmpsConfig& config) {
    VectorOperatorResult result;
    result.divergence.assign(particles.size(), 0.0);
    result.laplacian_x.assign(particles.size(), 0.0);
    result.laplacian_y.assign(particles.size(), 0.0);
    result.laplacian_z.assign(particles.size(), 0.0);
    result.divergence_error.assign(particles.size(), 0.0);
    result.laplacian_error.assign(particles.size(), 0.0);
    result.grad_xx.assign(particles.size(), 0.0);
    result.grad_xy.assign(particles.size(), 0.0);
    result.grad_xz.assign(particles.size(), 0.0);
    result.grad_yx.assign(particles.size(), 0.0);
    result.grad_yy.assign(particles.size(), 0.0);
    result.grad_yz.assign(particles.size(), 0.0);
    result.grad_zx.assign(particles.size(), 0.0);
    result.grad_zy.assign(particles.size(), 0.0);
    result.grad_zz.assign(particles.size(), 0.0);

    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!particles.isFluid(i) || matrices.particles[i].regular.status != lsmps::LsmpsMatrixStatus::Valid) {
            continue;
        }

        const auto& inverse = matrices.particles[i].regular.inverse_moment;
        const lsmps::LsmpsBasisVector dx = inverse * scalarRhs(particles, neighbors, i, ux, config);
        const lsmps::LsmpsBasisVector dy = inverse * scalarRhs(particles, neighbors, i, uy, config);
        const lsmps::LsmpsBasisVector dz = inverse * scalarRhs(particles, neighbors, i, uz, config);

        result.grad_xx[i] = dx[0] / support_radius;
        result.grad_xy[i] = dx[1] / support_radius;
        result.grad_xz[i] = dx[2] / support_radius;
        result.grad_yx[i] = dy[0] / support_radius;
        result.grad_yy[i] = dy[1] / support_radius;
        result.grad_yz[i] = dy[2] / support_radius;
        result.grad_zx[i] = dz[0] / support_radius;
        result.grad_zy[i] = dz[1] / support_radius;
        result.grad_zz[i] = dz[2] / support_radius;

        result.divergence[i] = dx[0] / support_radius + dy[1] / support_radius + dz[2] / support_radius;
        result.laplacian_x[i] = laplacianFromD(dx);
        result.laplacian_y[i] = laplacianFromD(dy);
        result.laplacian_z[i] = laplacianFromD(dz);

        const double analytical_lap_x = 0.0;
        const double analytical_lap_y = 2.0;
        const double analytical_lap_z = 4.0;
        result.divergence_error[i] = std::abs(result.divergence[i]);
        result.laplacian_error[i] = std::sqrt(
            (result.laplacian_x[i] - analytical_lap_x) * (result.laplacian_x[i] - analytical_lap_x) +
            (result.laplacian_y[i] - analytical_lap_y) * (result.laplacian_y[i] - analytical_lap_y) +
            (result.laplacian_z[i] - analytical_lap_z) * (result.laplacian_z[i] - analytical_lap_z));
    }
    return result;
}

void assertMaxBelow(
    const std::string& label,
    const std::vector<double>& values,
    const lsmps::ParticleSet& particles,
    double threshold) {
    double max_value = 0.0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            max_value = std::max(max_value, values[i]);
        }
    }
    if (max_value >= threshold) {
        std::cerr << label << " max error " << max_value << " exceeds threshold " << threshold << '\n';
    }
    assert(max_value < threshold);
}

void runHydrostaticCase() {
    lsmps::ParticleSet particles;
    particles.reserve(12000);
    addFluidBlock(particles, grid_size, grid_size, grid_size);
    addBoxWalls(particles, grid_size, grid_size, grid_size);

    std::vector<double>& pressure = particles.pressures();
    const double water_height = static_cast<double>(grid_size - 1) * spacing;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        pressure[i] = density * gravity_magnitude * (water_height - particles.positions()[i].y);
    }

    lsmps::NeighborSearch search(support_radius, {-spacing, -spacing, -spacing});
    const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
    search.updateNeighborCounts(particles, neighbors);
    const lsmps::LsmpsConfig config = validationConfig();
    const lsmps::LsmpsMatrixSet matrices =
        lsmps::buildLsmpsMatrices(particles, neighbors, support_radius, config);
    const ScalarOperatorResult result = computePressureOperators(particles, neighbors, matrices, pressure, config);

    assertMaxBelow("hydrostatic gradient", result.gradient_error, particles, 1.0e-8);
    assertMaxBelow("hydrostatic laplacian", result.laplacian_error, particles, 1.0e-8);

    const lsmps::VtkWriter writer;
    writer.writeParticles(
        "output/lsmps_hydrostatic_20x20x20.vtk",
        particles,
        {
            {"pressure_gradient_error", result.gradient_error},
            {"pressure_laplacian", result.laplacian},
            {"pressure_laplacian_error", result.laplacian_error},
            {"lsmps_pressure_status", matrixStatusField(particles, matrices, true)},
            {"lsmps_pressure_rank", matrixRankField(particles, matrices, true)},
            {"lsmps_pressure_condition_number", matrixConditionField(particles, matrices, true)},
            {"lsmps_wall_neighbor_count", wallCountField(particles, matrices)},
        },
        {
            {"pressure_gradient", result.gradient},
        });
}

void runPipeVelocityCase() {
    lsmps::ParticleSet particles;
    particles.reserve(12000);
    addFluidBlock(particles, grid_size, grid_size, grid_size);
    addBoxWalls(particles, grid_size, grid_size, grid_size);

    std::vector<double> ux(particles.size(), 0.0);
    std::vector<double> uy(particles.size(), 0.0);
    std::vector<double> uz(particles.size(), 0.0);
    std::vector<lsmps::Vector3> velocity(particles.size(), lsmps::Vector3{});
    for (std::size_t i = 0; i < particles.size(); ++i) {
        const lsmps::Vector3& x = particles.positions()[i];
        ux[i] = x.y * x.y - x.z * x.z;
        uy[i] = 2.0 * x.z * x.z - x.x * x.x;
        uz[i] = 3.0 * x.x * x.x - x.y * x.y;
        velocity[i] = {ux[i], uy[i], uz[i]};
        particles.velocities()[i] = velocity[i];
    }

    lsmps::NeighborSearch search(support_radius, {-spacing, -spacing, -spacing});
    const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
    search.updateNeighborCounts(particles, neighbors);
    const lsmps::LsmpsConfig config = validationConfig();
    const lsmps::LsmpsMatrixSet matrices =
        lsmps::buildLsmpsMatrices(particles, neighbors, support_radius, config);
    const VectorOperatorResult result = computeVelocityOperators(particles, neighbors, matrices, ux, uy, uz, config);

    assertMaxBelow("pipe divergence", result.divergence_error, particles, 1.0e-8);
    assertMaxBelow("pipe laplacian", result.laplacian_error, particles, 1.0e-8);

    std::vector<lsmps::Vector3> velocity_laplacian(particles.size(), lsmps::Vector3{});
    for (std::size_t i = 0; i < particles.size(); ++i) {
        velocity_laplacian[i] = {result.laplacian_x[i], result.laplacian_y[i], result.laplacian_z[i]};
    }

    const lsmps::VtkWriter writer;
    writer.writeParticles(
        "output/lsmps_pipe_20x20x20.vtk",
        particles,
        {
            {"velocity_magnitude", vectorMagnitude(velocity)},
            {"velocity_grad_xx", result.grad_xx},
            {"velocity_grad_xy", result.grad_xy},
            {"velocity_grad_xz", result.grad_xz},
            {"velocity_grad_yx", result.grad_yx},
            {"velocity_grad_yy", result.grad_yy},
            {"velocity_grad_yz", result.grad_yz},
            {"velocity_grad_zx", result.grad_zx},
            {"velocity_grad_zy", result.grad_zy},
            {"velocity_grad_zz", result.grad_zz},
            {"velocity_divergence", result.divergence},
            {"velocity_divergence_error", result.divergence_error},
            {"velocity_laplacian_x", result.laplacian_x},
            {"velocity_laplacian_y", result.laplacian_y},
            {"velocity_laplacian_z", result.laplacian_z},
            {"velocity_laplacian_error", result.laplacian_error},
            {"lsmps_regular_status", matrixStatusField(particles, matrices, false)},
            {"lsmps_regular_rank", matrixRankField(particles, matrices, false)},
            {"lsmps_regular_condition_number", matrixConditionField(particles, matrices, false)},
            {"lsmps_wall_neighbor_count", wallCountField(particles, matrices)},
        },
        {
            {"velocity_laplacian", velocity_laplacian},
        });
}

}  // namespace

int main() {
    runHydrostaticCase();
    runPipeVelocityCase();
    return 0;
}
