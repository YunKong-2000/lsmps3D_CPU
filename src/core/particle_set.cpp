#include "core/particle_set.hpp"

#include <stdexcept>

namespace lsmps {

std::size_t ParticleSet::size() const noexcept {
    return positions_.size();
}

bool ParticleSet::empty() const noexcept {
    return positions_.empty();
}

void ParticleSet::clear() {
    positions_.clear();
    velocities_.clear();
    pressures_.clear();
    densities_.clear();
    types_.clear();
    fluid_states_.clear();
    neighbor_counts_.clear();
}

void ParticleSet::reserve(std::size_t count) {
    positions_.reserve(count);
    velocities_.reserve(count);
    pressures_.reserve(count);
    densities_.reserve(count);
    types_.reserve(count);
    fluid_states_.reserve(count);
    neighbor_counts_.reserve(count);
}

std::size_t ParticleSet::addParticle(
    ParticleType type,
    const Vector3& position,
    const Vector3& velocity,
    double pressure,
    double density,
    FluidParticleState fluid_state) {
    const std::size_t index = size();
    positions_.push_back(position);
    velocities_.push_back(velocity);
    pressures_.push_back(pressure);
    densities_.push_back(density);
    types_.push_back(type);
    fluid_states_.push_back(fluid_state);
    neighbor_counts_.push_back(0);
    return index;
}

std::size_t ParticleSet::addFluidParticle(
    const Vector3& position,
    const Vector3& velocity,
    FluidParticleState fluid_state,
    double pressure,
    double density) {
    return addParticle(ParticleType::Fluid, position, velocity, pressure, density, fluid_state);
}

std::size_t ParticleSet::addWallParticle(
    const Vector3& position,
    const Vector3& velocity,
    double pressure,
    double density) {
    return addParticle(
        ParticleType::Wall,
        position,
        velocity,
        pressure,
        density,
        FluidParticleState::Internal);
}

const std::vector<Vector3>& ParticleSet::positions() const noexcept {
    return positions_;
}

std::vector<Vector3>& ParticleSet::positions() noexcept {
    return positions_;
}

const std::vector<Vector3>& ParticleSet::velocities() const noexcept {
    return velocities_;
}

std::vector<Vector3>& ParticleSet::velocities() noexcept {
    return velocities_;
}

const std::vector<double>& ParticleSet::pressures() const noexcept {
    return pressures_;
}

std::vector<double>& ParticleSet::pressures() noexcept {
    return pressures_;
}

const std::vector<double>& ParticleSet::densities() const noexcept {
    return densities_;
}

std::vector<double>& ParticleSet::densities() noexcept {
    return densities_;
}

const std::vector<ParticleType>& ParticleSet::types() const noexcept {
    return types_;
}

std::vector<ParticleType>& ParticleSet::types() noexcept {
    return types_;
}

const std::vector<FluidParticleState>& ParticleSet::fluidStates() const noexcept {
    return fluid_states_;
}

std::vector<FluidParticleState>& ParticleSet::fluidStates() noexcept {
    return fluid_states_;
}

const std::vector<std::size_t>& ParticleSet::neighborCounts() const noexcept {
    return neighbor_counts_;
}

std::vector<std::size_t>& ParticleSet::neighborCounts() noexcept {
    return neighbor_counts_;
}

bool ParticleSet::isFluid(std::size_t index) const {
    if (index >= size()) {
        throw std::out_of_range("ParticleSet::isFluid index out of range");
    }
    return types_[index] == ParticleType::Fluid;
}

bool ParticleSet::isWall(std::size_t index) const {
    if (index >= size()) {
        throw std::out_of_range("ParticleSet::isWall index out of range");
    }
    return types_[index] == ParticleType::Wall;
}

}  // namespace lsmps
