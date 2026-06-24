#include "core/particle_set.hpp"
#include "diagnostics/particle_distribution_diagnostics.hpp"
#include "neighbor/neighbor_search.hpp"

#include <cassert>
#include <cmath>
#include <cstddef>

namespace {

constexpr double spacing = 1.0;

std::size_t addGrid(lsmps::ParticleSet& particles, int nx, int ny, int nz) {
    std::size_t center = 0;
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const std::size_t index = particles.addFluidParticle(
                    {spacing * x, spacing * y, spacing * z},
                    {},
                    lsmps::FluidParticleState::Internal,
                    0.0,
                    1000.0);
                if (x == nx / 2 && y == ny / 2 && z == nz / 2) {
                    center = index;
                }
            }
        }
    }
    return center;
}

}  // namespace

int main() {
    lsmps::ParticleSet particles;
    const std::size_t center = addGrid(particles, 5, 5, 5);

    lsmps::NeighborSearch search(2.1 * spacing, {-spacing, -spacing, -spacing});
    const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
    search.updateNeighborCounts(particles, neighbors);

    lsmps::GeometryConfig geometry;
    geometry.particle_spacing = spacing;
    geometry.support_radius = 2.1 * spacing;

    const lsmps::ParticleDistributionDiagnostics diagnostics =
        lsmps::computeParticleDistributionDiagnostics(
            particles,
            neighbors,
            geometry,
            lsmps::LsmpsKernelType::Linear);

    assert(std::abs(diagnostics.nearest_fluid_distance[center] - spacing) < 1.0e-12);
    assert(diagnostics.max_fluid_neighbor_distance[center] > spacing);
    assert(diagnostics.mean_fluid_neighbor_distance[center] > spacing);
    assert(diagnostics.number_density[center] > 0.0);
    assert(diagnostics.number_density_ratio[center] > 0.99);
    assert(diagnostics.x_positive_neighbor_count[center] > 0.0);
    assert(diagnostics.x_negative_neighbor_count[center] > 0.0);
    assert(diagnostics.y_positive_neighbor_count[center] > 0.0);
    assert(diagnostics.y_negative_neighbor_count[center] > 0.0);
    assert(diagnostics.z_positive_neighbor_count[center] > 0.0);
    assert(diagnostics.z_negative_neighbor_count[center] > 0.0);
    assert(diagnostics.directional_coverage_score[center] > 0.0);
    assert(diagnostics.geometry_min_eigenvalue[center] > 0.0);
    assert(diagnostics.geometry_max_eigenvalue[center] > diagnostics.geometry_min_eigenvalue[center]);
    assert(diagnostics.geometry_condition_number[center] > 1.0);

    return 0;
}
