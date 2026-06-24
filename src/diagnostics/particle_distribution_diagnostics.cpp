#include "diagnostics/particle_distribution_diagnostics.hpp"

#include "lsmps/weight_function.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cstddef>
#include <limits>

namespace lsmps {
namespace {

constexpr double distance_epsilon = 1.0e-14;

double safeConditionNumber(double min_value, double max_value) {
    if (max_value <= distance_epsilon) {
        return 0.0;
    }
    if (min_value <= distance_epsilon) {
        return std::numeric_limits<double>::infinity();
    }
    return max_value / min_value;
}

}  // namespace

ParticleDistributionDiagnostics computeParticleDistributionDiagnostics(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    const GeometryConfig& geometry,
    LsmpsKernelType kernel_type) {
    const std::size_t count = particles.size();
    ParticleDistributionDiagnostics diagnostics;
    diagnostics.nearest_fluid_distance.assign(count, 0.0);
    diagnostics.max_fluid_neighbor_distance.assign(count, 0.0);
    diagnostics.mean_fluid_neighbor_distance.assign(count, 0.0);
    diagnostics.number_density.assign(count, 0.0);
    diagnostics.number_density_ratio.assign(count, 0.0);
    diagnostics.x_positive_neighbor_count.assign(count, 0.0);
    diagnostics.x_negative_neighbor_count.assign(count, 0.0);
    diagnostics.y_positive_neighbor_count.assign(count, 0.0);
    diagnostics.y_negative_neighbor_count.assign(count, 0.0);
    diagnostics.z_positive_neighbor_count.assign(count, 0.0);
    diagnostics.z_negative_neighbor_count.assign(count, 0.0);
    diagnostics.directional_coverage_score.assign(count, 0.0);
    diagnostics.geometry_min_eigenvalue.assign(count, 0.0);
    diagnostics.geometry_max_eigenvalue.assign(count, 0.0);
    diagnostics.geometry_condition_number.assign(count, 0.0);

    double reference_number_density = 0.0;
    const double support_radius = geometry.support_radius;

    for (std::size_t i = 0; i < count; ++i) {
        if (!particles.isFluid(i)) {
            continue;
        }

        const Vector3& position_i = particles.positions()[i];
        double nearest_distance = std::numeric_limits<double>::infinity();
        double max_distance = 0.0;
        double distance_sum = 0.0;
        std::size_t valid_fluid_neighbors = 0;
        Eigen::Matrix3d geometry_tensor = Eigen::Matrix3d::Zero();

        for (const std::size_t j : neighbors.fluid[i]) {
            const Vector3 offset = particles.positions()[j] - position_i;
            const double distance = norm(offset);
            if (distance <= distance_epsilon || distance > support_radius) {
                continue;
            }

            ++valid_fluid_neighbors;
            nearest_distance = std::min(nearest_distance, distance);
            max_distance = std::max(max_distance, distance);
            distance_sum += distance;

            const double weight = evaluateWeight(distance, support_radius, kernel_type);
            diagnostics.number_density[i] += weight;

            const Vector3 direction = offset / distance;
            if (direction.x > 0.0) {
                diagnostics.x_positive_neighbor_count[i] += 1.0;
            } else if (direction.x < 0.0) {
                diagnostics.x_negative_neighbor_count[i] += 1.0;
            }
            if (direction.y > 0.0) {
                diagnostics.y_positive_neighbor_count[i] += 1.0;
            } else if (direction.y < 0.0) {
                diagnostics.y_negative_neighbor_count[i] += 1.0;
            }
            if (direction.z > 0.0) {
                diagnostics.z_positive_neighbor_count[i] += 1.0;
            } else if (direction.z < 0.0) {
                diagnostics.z_negative_neighbor_count[i] += 1.0;
            }

            const Eigen::Vector3d e(direction.x, direction.y, direction.z);
            geometry_tensor.noalias() += weight * (e * e.transpose());
        }

        if (valid_fluid_neighbors > 0) {
            diagnostics.nearest_fluid_distance[i] = nearest_distance;
            diagnostics.max_fluid_neighbor_distance[i] = max_distance;
            diagnostics.mean_fluid_neighbor_distance[i] =
                distance_sum / static_cast<double>(valid_fluid_neighbors);
        }

        diagnostics.directional_coverage_score[i] = std::min({
            diagnostics.x_positive_neighbor_count[i],
            diagnostics.x_negative_neighbor_count[i],
            diagnostics.y_positive_neighbor_count[i],
            diagnostics.y_negative_neighbor_count[i],
            diagnostics.z_positive_neighbor_count[i],
            diagnostics.z_negative_neighbor_count[i],
        });

        if (valid_fluid_neighbors > 0) {
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(geometry_tensor);
            if (solver.info() == Eigen::Success) {
                diagnostics.geometry_min_eigenvalue[i] = solver.eigenvalues().minCoeff();
                diagnostics.geometry_max_eigenvalue[i] = solver.eigenvalues().maxCoeff();
                diagnostics.geometry_condition_number[i] = safeConditionNumber(
                    diagnostics.geometry_min_eigenvalue[i],
                    diagnostics.geometry_max_eigenvalue[i]);
            }
        }

        if (neighbors.wall[i].empty()) {
            reference_number_density = std::max(reference_number_density, diagnostics.number_density[i]);
        }
    }

    if (reference_number_density > distance_epsilon) {
        for (std::size_t i = 0; i < count; ++i) {
            if (particles.isFluid(i)) {
                diagnostics.number_density_ratio[i] = diagnostics.number_density[i] / reference_number_density;
            }
        }
    }

    return diagnostics;
}

}  // namespace lsmps
