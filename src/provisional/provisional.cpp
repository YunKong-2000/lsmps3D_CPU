#include "provisional/provisional.hpp"

#include "lsmps/lsmps_basis.hpp"
#include "lsmps/weight_function.hpp"

#include <cstddef>

namespace lsmps {
namespace {

std::vector<Vector3> normalizedExternalAcceleration(const ParticleSet& particles, const std::vector<Vector3>& external) {
    if (external.empty()) {
        return std::vector<Vector3>(particles.size(), Vector3{});
    }
    if (external.size() == particles.size()) {
        return external;
    }
    return std::vector<Vector3>(particles.size(), Vector3{});
}

LsmpsBasisVector scalarRhs(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    std::size_t i,
    const std::vector<double>& values,
    double support_radius,
    LsmpsKernelType kernel_type) {
    LsmpsBasisVector rhs = LsmpsBasisVector::Zero();
    const Vector3 position_i = particles.positions()[i];

    for (const std::size_t j : neighbors.fluid[i]) {
        const Vector3 offset = particles.positions()[j] - position_i;
        const double weight = evaluateWeight(norm(offset), support_radius, kernel_type);
        rhs.noalias() += weight * evaluateTypeABasis(offset, support_radius) * (values[j] - values[i]);
    }

    for (const std::size_t j : neighbors.wall[i]) {
        const Vector3 offset = particles.positions()[j] - position_i;
        const double weight = evaluateWeight(norm(offset), support_radius, kernel_type);
        rhs.noalias() += weight * evaluateTypeABasis(offset, support_radius) * (values[j] - values[i]);
    }

    return rhs;
}

double laplacianFromD(const LsmpsBasisVector& d, double support_radius) {
    return (d[3] + d[4] + d[5]) / (support_radius * support_radius);
}

}  // namespace

ProvisionalVelocityResult ProvisionalVelocityCalculator::compute(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    const LsmpsMatrixSet& matrices,
    const SimulationConfig& config,
    const std::vector<Vector3>& external_acceleration) const {
    ProvisionalVelocityResult result;
    result.computed = true;
    result.provisional_velocity = particles.velocities();
    result.viscous_acceleration.assign(particles.size(), {});
    result.body_acceleration.assign(particles.size(), {});
    result.velocity_delta.assign(particles.size(), {});
    result.status.assign(particles.size(), ProvisionalParticleStatus::NonFluidOrWall);

    std::vector<double> velocity_x(particles.size(), 0.0);
    std::vector<double> velocity_y(particles.size(), 0.0);
    std::vector<double> velocity_z(particles.size(), 0.0);
    for (std::size_t i = 0; i < particles.size(); ++i) {
        velocity_x[i] = particles.velocities()[i].x;
        velocity_y[i] = particles.velocities()[i].y;
        velocity_z[i] = particles.velocities()[i].z;
    }

    const std::vector<Vector3> external = normalizedExternalAcceleration(particles, external_acceleration);

    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!particles.isFluid(i) && !particles.isWall(i)) {
            continue;
        }

        result.body_acceleration[i] = config.physical.gravity + external[i];

        if (particles.isFluid(i)) {
            const LsmpsInverseMatrix& regular = matrices.particles[i].regular;
            if (regular.status != LsmpsMatrixStatus::Valid && regular.status != LsmpsMatrixStatus::IllConditioned) {
                result.status[i] = ProvisionalParticleStatus::MatrixUnavailable;
                result.velocity_delta[i] = config.time.dt * result.body_acceleration[i];
                result.provisional_velocity[i] = particles.velocities()[i] + result.velocity_delta[i];
                continue;
            }

            const LsmpsBasisVector dx = regular.inverse_moment *
                                        scalarRhs(
                                            particles,
                                            neighbors,
                                            i,
                                            velocity_x,
                                            config.geometry.support_radius,
                                            config.lsmps.kernel_type);
            const LsmpsBasisVector dy = regular.inverse_moment *
                                        scalarRhs(
                                            particles,
                                            neighbors,
                                            i,
                                            velocity_y,
                                            config.geometry.support_radius,
                                            config.lsmps.kernel_type);
            const LsmpsBasisVector dz = regular.inverse_moment *
                                        scalarRhs(
                                            particles,
                                            neighbors,
                                            i,
                                            velocity_z,
                                            config.geometry.support_radius,
                                            config.lsmps.kernel_type);

            result.viscous_acceleration[i] = {
                config.physical.viscosity * laplacianFromD(dx, config.geometry.support_radius),
                config.physical.viscosity * laplacianFromD(dy, config.geometry.support_radius),
                config.physical.viscosity * laplacianFromD(dz, config.geometry.support_radius),
            };
        }

        result.velocity_delta[i] =
            config.time.dt * (result.viscous_acceleration[i] + result.body_acceleration[i]);
        result.provisional_velocity[i] = particles.velocities()[i] + result.velocity_delta[i];
        result.status[i] = ProvisionalParticleStatus::Updated;
    }

    return result;
}

}  // namespace lsmps
