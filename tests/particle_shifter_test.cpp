#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "neighbor/neighbor_search.hpp"
#include "particle_shifting/particle_shifter.hpp"

#include <cassert>
#include <cmath>

namespace {

double distance(const lsmps::Vector3& a, const lsmps::Vector3& b) {
    return lsmps::norm(a - b);
}

}  // namespace

int main() {
    lsmps::SimulationConfig config;
    config.geometry.particle_spacing = 1.0;
    config.geometry.support_radius = 2.0;
    config.particle_shifting.enabled = true;
    config.particle_shifting.max_displacement_factor = 0.1;
    config.particle_shifting.min_distance_factor = 0.7;

    lsmps::ParticleSet particles;
    const std::size_t a = particles.addFluidParticle(
        {0.0, 0.0, 0.0},
        {},
        lsmps::FluidParticleState::Internal);
    const std::size_t b = particles.addFluidParticle(
        {0.5, 0.0, 0.0},
        {},
        lsmps::FluidParticleState::Internal);
    const std::size_t c = particles.addFluidParticle(
        {0.0, 1.5, 0.0},
        {},
        lsmps::FluidParticleState::FreeSurface);

    lsmps::NeighborSearch search(config.geometry.support_radius, {-1.0, -1.0, -1.0});
    const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
    const double initial_distance = distance(particles.positions()[a], particles.positions()[b]);

    const lsmps::ParticleShifter shifter;
    const lsmps::ParticleShiftResult result = shifter.compute(particles, neighbors, config);
    assert(result.applied);
    assert(result.repulsion_active[a] == 1.0);
    assert(result.repulsion_active[b] == 1.0);
    assert(result.repulsion_active[c] == 0.0);
    assert(result.magnitude[a] > 0.0);
    assert(result.magnitude[a] <= 0.1 + 1.0e-12);
    assert(result.magnitude[c] == 0.0);

    shifter.apply(particles, result);
    const double shifted_distance = distance(particles.positions()[a], particles.positions()[b]);
    assert(shifted_distance > initial_distance);
    assert(particles.positions()[c].x == 0.0);
    assert(particles.positions()[c].y == 1.5);
    assert(particles.positions()[c].z == 0.0);

    config.particle_shifting.enabled = false;
    const lsmps::ParticleShiftResult disabled = shifter.compute(particles, neighbors, config);
    assert(!disabled.applied);
    for (double magnitude : disabled.magnitude) {
        assert(magnitude == 0.0);
    }

    return 0;
}
