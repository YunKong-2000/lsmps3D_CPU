#pragma once

#include "core/particle_types.hpp"
#include "core/vector3.hpp"

#include <cstddef>
#include <vector>

namespace lsmps {

class ParticleSet {
public:
    std::size_t size() const noexcept;
    bool empty() const noexcept;
    void clear();
    void reserve(std::size_t count);

    std::size_t addParticle(
        ParticleType type,
        const Vector3& position,
        const Vector3& velocity = {},
        double pressure = 0.0,
        double density = 0.0,
        FluidParticleState fluid_state = FluidParticleState::Internal,
        const Vector3& wall_normal = {});

    std::size_t addFluidParticle(
        const Vector3& position,
        const Vector3& velocity = {},
        FluidParticleState fluid_state = FluidParticleState::Internal,
        double pressure = 0.0,
        double density = 0.0);

    std::size_t addWallParticle(
        const Vector3& position,
        const Vector3& velocity = {},
        double pressure = 0.0,
        double density = 0.0,
        const Vector3& wall_normal = {});

    const std::vector<Vector3>& positions() const noexcept;
    std::vector<Vector3>& positions() noexcept;

    const std::vector<Vector3>& velocities() const noexcept;
    std::vector<Vector3>& velocities() noexcept;

    const std::vector<double>& pressures() const noexcept;
    std::vector<double>& pressures() noexcept;

    const std::vector<double>& densities() const noexcept;
    std::vector<double>& densities() noexcept;

    const std::vector<ParticleType>& types() const noexcept;
    std::vector<ParticleType>& types() noexcept;

    const std::vector<FluidParticleState>& fluidStates() const noexcept;
    std::vector<FluidParticleState>& fluidStates() noexcept;

    const std::vector<Vector3>& wallNormals() const noexcept;
    std::vector<Vector3>& wallNormals() noexcept;

    const std::vector<std::size_t>& neighborCounts() const noexcept;
    std::vector<std::size_t>& neighborCounts() noexcept;

    const std::vector<std::size_t>& fluidNeighborCounts() const noexcept;
    std::vector<std::size_t>& fluidNeighborCounts() noexcept;

    const std::vector<std::size_t>& wallNeighborCounts() const noexcept;
    std::vector<std::size_t>& wallNeighborCounts() noexcept;

    bool isFluid(std::size_t index) const;
    bool isWall(std::size_t index) const;

private:
    std::vector<Vector3> positions_;
    std::vector<Vector3> velocities_;
    std::vector<double> pressures_;
    std::vector<double> densities_;
    std::vector<ParticleType> types_;
    std::vector<FluidParticleState> fluid_states_;
    std::vector<Vector3> wall_normals_;
    std::vector<std::size_t> neighbor_counts_;
    std::vector<std::size_t> fluid_neighbor_counts_;
    std::vector<std::size_t> wall_neighbor_counts_;
};

}  // namespace lsmps
