#include "core/particle_set.hpp"
#include "core/simulation_state.hpp"
#include "core/vector3.hpp"

#include <cassert>
#include <cstddef>

int main() {
    lsmps::ParticleSet particles;
    particles.reserve(2);

    const std::size_t fluid = particles.addFluidParticle(
        {1.0, 2.0, 3.0},
        {0.1, 0.2, 0.3},
        lsmps::FluidParticleState::FreeSurface,
        100.0,
        1000.0);

    const std::size_t wall = particles.addWallParticle({0.0, 0.0, 0.0});

    assert(fluid == 0);
    assert(wall == 1);
    assert(particles.size() == 2);
    assert(particles.isFluid(fluid));
    assert(particles.isWall(wall));
    assert(particles.types()[fluid] == lsmps::ParticleType::Fluid);
    assert(particles.types()[wall] == lsmps::ParticleType::Wall);
    assert(particles.fluidStates()[fluid] == lsmps::FluidParticleState::FreeSurface);
    assert(particles.positions()[fluid].x == 1.0);
    assert(particles.velocities()[fluid].z == 0.3);
    assert(particles.pressures()[fluid] == 100.0);
    assert(particles.densities()[fluid] == 1000.0);
    assert(particles.neighborCounts()[fluid] == 0);
    assert(particles.fluidNeighborCounts()[fluid] == 0);
    assert(particles.wallNeighborCounts()[fluid] == 0);

    particles.pressures()[fluid] = 200.0;
    particles.neighborCounts()[fluid] = 12;
    particles.fluidNeighborCounts()[fluid] = 9;
    particles.wallNeighborCounts()[fluid] = 3;
    particles.fluidStates()[fluid] = lsmps::FluidParticleState::NearFreeSurface;

    assert(particles.pressures()[fluid] == 200.0);
    assert(particles.neighborCounts()[fluid] == 12);
    assert(particles.fluidNeighborCounts()[fluid] == 9);
    assert(particles.wallNeighborCounts()[fluid] == 3);
    assert(particles.fluidStates()[fluid] == lsmps::FluidParticleState::NearFreeSurface);

    const lsmps::Vector3 sum = particles.positions()[fluid] + lsmps::Vector3{1.0, 1.0, 1.0};
    assert(sum.x == 2.0);
    assert(sum.y == 3.0);
    assert(sum.z == 4.0);
    assert(lsmps::dot(sum, lsmps::Vector3{1.0, 0.0, 0.0}) == 2.0);

    lsmps::SimulationState state;
    state.current_step = 3;
    state.current_time = 0.15;
    state.max_neighbor_count = particles.neighborCounts()[fluid];

    assert(state.current_step == 3);
    assert(state.current_time == 0.15);
    assert(state.max_neighbor_count == 12);

    particles.clear();
    assert(particles.empty());

    return 0;
}
