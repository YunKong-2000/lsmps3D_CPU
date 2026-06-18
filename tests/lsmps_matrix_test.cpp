#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "core/vector3.hpp"
#include "io/vtk_writer.hpp"
#include "lsmps/lsmps_basis.hpp"
#include "lsmps/lsmps_matrices.hpp"
#include "lsmps/weight_function.hpp"
#include "neighbor/neighbor_search.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace {

constexpr double spacing = 1.0;
constexpr double support_radius = 3.1 * spacing;
constexpr double tolerance = 1.0e-8;

std::size_t gridIndex(int x, int y, int z, int nx, int ny) {
    return static_cast<std::size_t>((z * ny + y) * nx + x);
}

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

double polynomial(const lsmps::Vector3& position) {
    return position.x * position.x + 2.0 * position.y * position.y + 3.0 * position.z * position.z +
           position.x * position.y - 2.0 * position.x * position.z + position.y * position.z;
}

lsmps::LsmpsBasisVector analyticalD(const lsmps::Vector3& position) {
    lsmps::LsmpsBasisVector expected;
    expected << support_radius * (2.0 * position.x + position.y - 2.0 * position.z),
        support_radius * (4.0 * position.y + position.x + position.z),
        support_radius * (6.0 * position.z - 2.0 * position.x + position.y),
        support_radius * support_radius * 2.0,
        support_radius * support_radius * 4.0,
        support_radius * support_radius * 6.0,
        support_radius * support_radius,
        -2.0 * support_radius * support_radius,
        support_radius * support_radius;
    return expected;
}

std::vector<double> statusField(
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

std::vector<double> rankField(
    const lsmps::ParticleSet& particles,
    const lsmps::LsmpsMatrixSet& matrices,
    bool pressure) {
    std::vector<double> values(particles.size(), 0.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        values[i] = static_cast<double>(pressure ? matrices.particles[i].pressure_neumann.rank : matrices.particles[i].regular.rank);
    }
    return values;
}

std::vector<double> conditionField(
    const lsmps::ParticleSet& particles,
    const lsmps::LsmpsMatrixSet& matrices,
    bool pressure) {
    std::vector<double> values(particles.size(), 0.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        values[i] = pressure ? matrices.particles[i].pressure_neumann.condition_number : matrices.particles[i].regular.condition_number;
    }
    return values;
}

std::vector<double> countField(
    const lsmps::ParticleSet& particles,
    const lsmps::LsmpsMatrixSet& matrices,
    bool wall) {
    std::vector<double> values(particles.size(), 0.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        values[i] = static_cast<double>(
            wall ? matrices.particles[i].regular.wall_neighbor_count : matrices.particles[i].regular.fluid_neighbor_count);
    }
    return values;
}

void writeDiagnostics(
    const std::string& path,
    const lsmps::ParticleSet& particles,
    const lsmps::LsmpsMatrixSet& matrices) {
    const lsmps::VtkWriter writer;
    writer.writeParticles(
        path,
        particles,
        {
            {"lsmps_regular_status", statusField(particles, matrices, false)},
            {"lsmps_regular_rank", rankField(particles, matrices, false)},
            {"lsmps_regular_condition_number", conditionField(particles, matrices, false)},
            {"lsmps_pressure_status", statusField(particles, matrices, true)},
            {"lsmps_pressure_rank", rankField(particles, matrices, true)},
            {"lsmps_pressure_condition_number", conditionField(particles, matrices, true)},
            {"lsmps_fluid_neighbor_count", countField(particles, matrices, false)},
            {"lsmps_wall_neighbor_count", countField(particles, matrices, true)},
        });
}

}  // namespace

int main() {
    {
        const lsmps::LsmpsBasisVector p = lsmps::evaluateTypeABasis({1.0, 2.0, 3.0}, 2.0);
        assert(std::abs(p[0] - 0.5) < tolerance);
        assert(std::abs(p[3] - 0.125) < tolerance);
        assert(std::abs(p[6] - 0.5) < tolerance);

        const lsmps::LsmpsBasisVector q =
            lsmps::evaluateTypeANeumannBasis({1.0, 2.0, 3.0}, {1.0, 0.0, 0.0}, 2.0);
        assert(std::abs(q[0] - 1.0) < tolerance);
        assert(std::abs(q[3] - 0.5) < tolerance);
        assert(std::abs(q[6] - 1.0) < tolerance);
        assert(std::abs(q[7] - 1.5) < tolerance);
    }

    lsmps::ParticleSet particles;
    constexpr int nx = 7;
    constexpr int ny = 7;
    constexpr int nz = 7;
    particles.reserve(1000);
    addFluidBlock(particles, nx, ny, nz);
    addBoxWalls(particles, nx, ny, nz);

    lsmps::NeighborSearch search(support_radius, {-1.0, -1.0, -1.0});
    const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
    search.updateNeighborCounts(particles, neighbors);

    lsmps::LsmpsConfig config;
    config.min_neighbors = 9;
    config.eigenvalue_tolerance = 1.0e-12;
    config.condition_number_warning = 1.0e10;
    config.condition_number_failure = 1.0e14;

    const lsmps::LsmpsMatrixSet matrices =
        lsmps::buildLsmpsMatrices(particles, neighbors, support_radius, config);
    writeDiagnostics("output/lsmps_matrix_diagnostics.vtk", particles, matrices);

    const std::size_t center = gridIndex(3, 3, 3, nx, ny);
    const auto& center_matrices = matrices.particles[center];
    assert(center_matrices.regular.status == lsmps::LsmpsMatrixStatus::Valid);
    assert(center_matrices.pressure_neumann.status == lsmps::LsmpsMatrixStatus::Valid);
    assert(center_matrices.regular.wall_neighbor_count == 0);
    assert((center_matrices.regular.inverse_moment - center_matrices.pressure_neumann.inverse_moment).norm() < tolerance);

    lsmps::LsmpsBasisVector rhs = lsmps::LsmpsBasisVector::Zero();
    for (const std::size_t j : neighbors.fluid[center]) {
        const lsmps::Vector3 offset = particles.positions()[j] - particles.positions()[center];
        const double distance = lsmps::norm(offset);
        const double weight = lsmps::evaluateWeight(distance, support_radius, config.kernel_type);
        rhs.noalias() += weight * lsmps::evaluateTypeABasis(offset, support_radius) *
                         (polynomial(particles.positions()[j]) - polynomial(particles.positions()[center]));
    }
    for (const std::size_t j : neighbors.wall[center]) {
        const lsmps::Vector3 offset = particles.positions()[j] - particles.positions()[center];
        const double distance = lsmps::norm(offset);
        const double weight = lsmps::evaluateWeight(distance, support_radius, config.kernel_type);
        rhs.noalias() += weight * lsmps::evaluateTypeABasis(offset, support_radius) *
                         (polynomial(particles.positions()[j]) - polynomial(particles.positions()[center]));
    }
    const lsmps::LsmpsBasisVector reconstructed = center_matrices.regular.inverse_moment * rhs;
    assert((reconstructed - analyticalD(particles.positions()[center])).norm() < 1.0e-7);

    const std::size_t near_wall = gridIndex(0, 3, 3, nx, ny);
    const auto& near_wall_matrices = matrices.particles[near_wall];
    assert(near_wall_matrices.regular.status == lsmps::LsmpsMatrixStatus::Valid);
    assert(near_wall_matrices.pressure_neumann.status == lsmps::LsmpsMatrixStatus::Valid);
    assert(near_wall_matrices.regular.wall_neighbor_count > 0);
    assert((near_wall_matrices.regular.inverse_moment - near_wall_matrices.pressure_neumann.inverse_moment).norm() > 1.0e-6);

    lsmps::LsmpsBasisVector pressure_rhs = lsmps::LsmpsBasisVector::Zero();
    const lsmps::Vector3 constant_gradient = {2.0, -3.0, 4.0};
    for (const std::size_t j : neighbors.fluid[near_wall]) {
        const lsmps::Vector3 offset = particles.positions()[j] - particles.positions()[near_wall];
        const double distance = lsmps::norm(offset);
        const double weight = lsmps::evaluateWeight(distance, support_radius, config.kernel_type);
        const double pressure_difference = lsmps::dot(constant_gradient, offset);
        pressure_rhs.noalias() += weight * lsmps::evaluateTypeABasis(offset, support_radius) * pressure_difference;
    }
    for (const std::size_t j : neighbors.wall[near_wall]) {
        const lsmps::Vector3 offset = particles.positions()[j] - particles.positions()[near_wall];
        const double distance = lsmps::norm(offset);
        const double weight = lsmps::evaluateWeight(distance, support_radius, config.kernel_type);
        const double normal_derivative = lsmps::dot(constant_gradient, particles.wallNormals()[j]);
        pressure_rhs.noalias() += weight *
                                  lsmps::evaluateTypeANeumannBasis(offset, particles.wallNormals()[j], support_radius) *
                                  support_radius * normal_derivative;
    }
    const lsmps::LsmpsBasisVector pressure_d =
        near_wall_matrices.pressure_neumann.inverse_moment * pressure_rhs;
    assert(std::abs(pressure_d[0] / support_radius - constant_gradient.x) < 1.0e-8);
    assert(std::abs(pressure_d[1] / support_radius - constant_gradient.y) < 1.0e-8);
    assert(std::abs(pressure_d[2] / support_radius - constant_gradient.z) < 1.0e-8);

    return 0;
}
