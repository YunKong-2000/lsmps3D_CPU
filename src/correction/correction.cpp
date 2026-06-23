#include "correction/correction.hpp"

#include "lsmps/lsmps_basis.hpp"
#include "lsmps/weight_function.hpp"

#include <cstddef>
#include <stdexcept>

namespace lsmps {
namespace {

LsmpsBasisVector pressureRhs(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    std::size_t i,
    const std::vector<double>& pressure,
    const SimulationConfig& config) {
    LsmpsBasisVector rhs = LsmpsBasisVector::Zero();
    const Vector3 position_i = particles.positions()[i];
    const double re = config.geometry.support_radius;

    for (const std::size_t j : neighbors.fluid[i]) {
        const Vector3 offset = particles.positions()[j] - position_i;
        const double weight = evaluateWeight(norm(offset), re, config.lsmps.kernel_type);
        rhs.noalias() += weight * evaluateTypeABasis(offset, re) * (pressure[j] - pressure[i]);
    }

    for (const std::size_t j : neighbors.wall[i]) {
        const Vector3 offset = particles.positions()[j] - position_i;
        const double weight = evaluateWeight(norm(offset), re, config.lsmps.kernel_type);
        const double neumann = config.physical.density * dot(config.physical.gravity, particles.wallNormals()[j]);
        rhs.noalias() += weight * evaluateTypeANeumannBasis(offset, particles.wallNormals()[j], re) * re * neumann;
    }

    return rhs;
}

Vector3 gradientFromD(const LsmpsBasisVector& d, double support_radius) {
    return {d[0] / support_radius, d[1] / support_radius, d[2] / support_radius};
}

}  // namespace

CorrectionResult PressureCorrectionApplier::apply(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    const LsmpsMatrixSet& matrices,
    const ProvisionalVelocityResult& provisional,
    const PressurePoissonResult& pressure,
    const SimulationConfig& config) const {
    if (!provisional.computed || provisional.provisional_velocity.size() != particles.size()) {
        throw std::runtime_error("PressureCorrectionApplier requires computed provisional velocities");
    }
    if (!pressure.solved || pressure.pressure.size() != particles.size()) {
        throw std::runtime_error("PressureCorrectionApplier requires solved pressure values");
    }

    CorrectionResult result;
    result.applied = true;
    result.pressure_gradient.assign(particles.size(), {});
    result.velocity_correction.assign(particles.size(), {});
    result.next_velocity = particles.velocities();
    result.displacement.assign(particles.size(), {});
    result.next_position = particles.positions();
    result.status.assign(particles.size(), CorrectionParticleStatus::NonFluid);

    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!particles.isFluid(i)) {
            continue;
        }

        const LsmpsInverseMatrix& pressure_matrix = matrices.particles[i].pressure_neumann;
        if (pressure_matrix.status != LsmpsMatrixStatus::Valid &&
            pressure_matrix.status != LsmpsMatrixStatus::IllConditioned) {
            result.status[i] = CorrectionParticleStatus::MatrixUnavailable;
            result.next_velocity[i] = provisional.provisional_velocity[i];
            result.displacement[i] =
                0.5 * config.time.dt * (particles.velocities()[i] + result.next_velocity[i]);
            result.next_position[i] = particles.positions()[i] + result.displacement[i];
            continue;
        }

        const LsmpsBasisVector d =
            pressure_matrix.inverse_moment * pressureRhs(particles, neighbors, i, pressure.pressure, config);
        result.pressure_gradient[i] = gradientFromD(d, config.geometry.support_radius);
        result.velocity_correction[i] =
            -(config.time.dt / config.physical.density) * result.pressure_gradient[i];
        result.next_velocity[i] = provisional.provisional_velocity[i] + result.velocity_correction[i];
        result.displacement[i] =
            0.5 * config.time.dt * (particles.velocities()[i] + result.next_velocity[i]);
        result.next_position[i] = particles.positions()[i] + result.displacement[i];
        result.status[i] = CorrectionParticleStatus::Updated;
    }

    return result;
}

}  // namespace lsmps
