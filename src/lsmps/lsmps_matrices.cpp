#include "lsmps/lsmps_matrices.hpp"

#include "lsmps/weight_function.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace lsmps {
namespace {

constexpr double distance_epsilon = 1.0e-14;

struct RawBuildResult {
    LsmpsMomentMatrix raw = LsmpsMomentMatrix::Zero();
    std::size_t fluid_count = 0;
    std::size_t wall_count = 0;
};

void accumulateOuter(
    RawBuildResult& result,
    const LsmpsBasisVector& basis,
    double weight,
    bool is_wall) {
    result.raw.noalias() += weight * (basis * basis.transpose());
    if (is_wall) {
        ++result.wall_count;
    } else {
        ++result.fluid_count;
    }
}

bool validNeighbor(
    const Vector3& offset,
    double support_radius,
    double& distance,
    double& weight,
    LsmpsKernelType kernel_type) {
    distance = norm(offset);
    if (distance <= distance_epsilon || distance > support_radius) {
        return false;
    }

    weight = evaluateWeight(distance, support_radius, kernel_type);
    return weight > 0.0;
}

RawBuildResult buildRegularRaw(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    std::size_t i,
    double support_radius,
    LsmpsKernelType kernel_type) {
    RawBuildResult result;
    const Vector3& position_i = particles.positions()[i];

    for (const std::size_t j : neighbors.fluid[i]) {
        double distance = 0.0;
        double weight = 0.0;
        const Vector3 offset = particles.positions()[j] - position_i;
        if (!validNeighbor(offset, support_radius, distance, weight, kernel_type)) {
            continue;
        }
        (void)distance;
        accumulateOuter(result, evaluateTypeABasis(offset, support_radius), weight, false);
    }

    for (const std::size_t j : neighbors.wall[i]) {
        double distance = 0.0;
        double weight = 0.0;
        const Vector3 offset = particles.positions()[j] - position_i;
        if (!validNeighbor(offset, support_radius, distance, weight, kernel_type)) {
            continue;
        }
        (void)distance;
        accumulateOuter(result, evaluateTypeABasis(offset, support_radius), weight, true);
    }

    return result;
}

RawBuildResult buildPressureRaw(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    std::size_t i,
    double support_radius,
    LsmpsKernelType kernel_type) {
    RawBuildResult result;
    const Vector3& position_i = particles.positions()[i];

    for (const std::size_t j : neighbors.fluid[i]) {
        double distance = 0.0;
        double weight = 0.0;
        const Vector3 offset = particles.positions()[j] - position_i;
        if (!validNeighbor(offset, support_radius, distance, weight, kernel_type)) {
            continue;
        }
        (void)distance;
        accumulateOuter(result, evaluateTypeABasis(offset, support_radius), weight, false);
    }

    for (const std::size_t j : neighbors.wall[i]) {
        double distance = 0.0;
        double weight = 0.0;
        const Vector3 offset = particles.positions()[j] - position_i;
        if (!validNeighbor(offset, support_radius, distance, weight, kernel_type)) {
            continue;
        }
        (void)distance;
        accumulateOuter(
            result,
            evaluateTypeANeumannBasis(offset, particles.wallNormals()[j], support_radius),
            weight,
            true);
    }

    return result;
}

LsmpsInverseMatrix invertAndDiagnose(const RawBuildResult& raw_result, const LsmpsConfig& config) {
    LsmpsInverseMatrix output;
    output.fluid_neighbor_count = raw_result.fluid_count;
    output.wall_neighbor_count = raw_result.wall_count;
    output.total_neighbor_count = raw_result.fluid_count + raw_result.wall_count;

    if (output.total_neighbor_count < config.min_neighbors) {
        output.status = LsmpsMatrixStatus::NotEnoughNeighbors;
        return output;
    }

    Eigen::SelfAdjointEigenSolver<LsmpsMomentMatrix> eigen_solver(raw_result.raw);
    if (eigen_solver.info() != Eigen::Success) {
        output.status = LsmpsMatrixStatus::InversionFailed;
        return output;
    }

    const auto eigenvalues = eigen_solver.eigenvalues();
    output.min_eigenvalue = eigenvalues.minCoeff();
    output.max_eigenvalue = eigenvalues.maxCoeff();
    output.rank = 0;
    for (int k = 0; k < eigenvalues.size(); ++k) {
        if (eigenvalues[k] > config.eigenvalue_tolerance) {
            ++output.rank;
        }
    }

    if (output.max_eigenvalue <= config.eigenvalue_tolerance) {
        output.condition_number = std::numeric_limits<double>::infinity();
    } else {
        const double denominator = std::max(output.min_eigenvalue, config.eigenvalue_tolerance);
        output.condition_number = output.max_eigenvalue / denominator;
    }

    if (output.rank < lsmps_basis_size || output.min_eigenvalue <= config.eigenvalue_tolerance) {
        output.status = LsmpsMatrixStatus::RankDeficient;
        return output;
    }

    Eigen::FullPivLU<LsmpsMomentMatrix> lu(raw_result.raw);
    if (!lu.isInvertible()) {
        output.status = LsmpsMatrixStatus::InversionFailed;
        return output;
    }

    output.inverse_moment = lu.inverse();
    if (!output.inverse_moment.allFinite()) {
        output.inverse_moment.setZero();
        output.status = LsmpsMatrixStatus::InversionFailed;
        return output;
    }

    if (output.condition_number > config.condition_number_failure) {
        output.status = LsmpsMatrixStatus::InversionFailed;
    } else if (output.condition_number > config.condition_number_warning) {
        output.status = LsmpsMatrixStatus::IllConditioned;
    } else {
        output.status = LsmpsMatrixStatus::Valid;
    }
    return output;
}

}  // namespace

LsmpsMatrixSet buildLsmpsMatrices(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    double support_radius,
    const LsmpsConfig& config) {
    LsmpsMatrixSet matrices;
    matrices.particles.resize(particles.size());

    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!particles.isFluid(i)) {
            continue;
        }

        const RawBuildResult regular_raw =
            buildRegularRaw(particles, neighbors, i, support_radius, config.kernel_type);
        matrices.particles[i].regular = invertAndDiagnose(regular_raw, config);

        const RawBuildResult pressure_raw =
            buildPressureRaw(particles, neighbors, i, support_radius, config.kernel_type);
        matrices.particles[i].pressure_neumann = invertAndDiagnose(pressure_raw, config);
    }

    return matrices;
}

}  // namespace lsmps
