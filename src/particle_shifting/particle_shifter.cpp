#include "particle_shifting/particle_shifter.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace lsmps {
namespace {

bool canShiftFluidParticle(const ParticleSet& particles, std::size_t index) {
    if (!particles.isFluid(index)) {
        return false;
    }
    const FluidParticleState state = particles.fluidStates()[index];
    return state == FluidParticleState::Internal || state == FluidParticleState::NearFreeSurface;
}

Vector3 limitedVector(const Vector3& value, double limit, bool& limited) {
    const double length = norm(value);
    if (limit <= 0.0 || length <= limit) {
        return value;
    }
    limited = true;
    return (limit / length) * value;
}

}  // namespace

ParticleShiftResult ParticleShifter::compute(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors,
    const SimulationConfig& config) const {
    ParticleShiftResult result;
    result.applied = config.particle_shifting.enabled;
    result.displacement.assign(particles.size(), {});
    result.magnitude.assign(particles.size(), 0.0);
    result.limited.assign(particles.size(), 0.0);
    result.repulsion_active.assign(particles.size(), 0.0);

    if (!config.particle_shifting.enabled) {
        return result;
    }

    const double h = config.geometry.particle_spacing;
    const double min_distance = config.particle_shifting.min_distance_factor * h;
    const double max_displacement = config.particle_shifting.max_displacement_factor * h;
    const double epsilon = std::max(1.0e-12 * h, std::numeric_limits<double>::epsilon());

    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!canShiftFluidParticle(particles, i)) {
            continue;
        }

        Vector3 shift;
        for (const std::size_t j : neighbors.fluid[i]) {
            const Vector3 offset = particles.positions()[i] - particles.positions()[j];
            const double distance = norm(offset);
            if (distance <= epsilon || distance >= min_distance) {
                continue;
            }

            const double strength = (min_distance - distance) / min_distance;
            shift += strength * (offset / (distance + epsilon));
            result.repulsion_active[i] = 1.0;
        }

        if (result.repulsion_active[i] == 0.0) {
            continue;
        }

        shift *= max_displacement;
        bool was_limited = false;
        shift = limitedVector(shift, max_displacement, was_limited);
        result.displacement[i] = shift;
        result.magnitude[i] = norm(shift);
        result.limited[i] = was_limited ? 1.0 : 0.0;
    }

    return result;
}

void ParticleShifter::apply(ParticleSet& particles, const ParticleShiftResult& result) const {
    if (!result.applied) {
        return;
    }
    if (result.displacement.size() != particles.size()) {
        throw std::runtime_error("ParticleShiftResult displacement size does not match particle count");
    }

    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            particles.positions()[i] += result.displacement[i];
        }
    }
}

}  // namespace lsmps
