#include "core/particle_set.hpp"
#include "core/vector3.hpp"
#include "neighbor/cell_linked_list.hpp"
#include "neighbor/neighbor_search.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace {

lsmps::NeighborList buildBruteForceNeighbors(const lsmps::ParticleSet& particles, double support_radius) {
    lsmps::NeighborList neighbors(particles.size());
    const double support_radius_squared = support_radius * support_radius;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        for (std::size_t j = 0; j < particles.size(); ++j) {
            if (i == j) {
                continue;
            }
            const lsmps::Vector3 offset = particles.positions()[j] - particles.positions()[i];
            if (lsmps::normSquared(offset) <= support_radius_squared) {
                neighbors[i].push_back(j);
            }
        }
        std::sort(neighbors[i].begin(), neighbors[i].end());
    }
    return neighbors;
}

bool contains(const std::vector<std::size_t>& values, std::size_t target) {
    return std::find(values.begin(), values.end(), target) != values.end();
}

std::size_t gridIndex(int x, int y, int z, int nx, int ny) {
    return static_cast<std::size_t>((z * ny + y) * nx + x);
}

}  // namespace

int main() {
    lsmps::ParticleSet grid_particles;
    for (int z = 0; z < 3; ++z) {
        for (int y = 0; y < 3; ++y) {
            for (int x = 0; x < 3; ++x) {
                grid_particles.addFluidParticle({static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)});
            }
        }
    }

    lsmps::NeighborSearch grid_search(1.01, {0.0, 0.0, 0.0});
    const lsmps::NeighborList grid_neighbors = grid_search.buildNeighborList(grid_particles);
    const lsmps::NeighborList brute_force_neighbors = buildBruteForceNeighbors(grid_particles, 1.01);
    assert(grid_neighbors == brute_force_neighbors);

    const std::size_t corner = 0;
    const std::size_t center = 13;
    assert(grid_neighbors[corner].size() == 3);
    assert(grid_neighbors[center].size() == 6);
    assert(contains(grid_neighbors[center], 4));
    assert(contains(grid_neighbors[center], 10));
    assert(contains(grid_neighbors[center], 12));
    assert(contains(grid_neighbors[center], 14));
    assert(contains(grid_neighbors[center], 16));
    assert(contains(grid_neighbors[center], 22));

    grid_search.updateNeighborCounts(grid_particles, grid_neighbors);
    assert(grid_particles.neighborCounts()[corner] == 3);
    assert(grid_particles.fluidNeighborCounts()[corner] == 3);
    assert(grid_particles.wallNeighborCounts()[corner] == 0);
    assert(grid_particles.neighborCounts()[center] == 6);
    assert(grid_particles.fluidNeighborCounts()[center] == 6);
    assert(grid_particles.wallNeighborCounts()[center] == 0);
    assert(grid_search.cells().cellCount() == 8);

    lsmps::ParticleSet mixed_particles;
    const std::size_t fluid = mixed_particles.addFluidParticle({0.0, 0.0, 0.0});
    const std::size_t wall = mixed_particles.addWallParticle({0.5, 0.0, 0.0});
    const std::size_t far_fluid = mixed_particles.addFluidParticle({2.0, 0.0, 0.0});

    lsmps::NeighborSearch mixed_search(0.75, {0.0, 0.0, 0.0});
    const lsmps::TypedNeighborList mixed_typed_neighbors = mixed_search.buildTypedNeighborList(mixed_particles);
    const lsmps::NeighborList mixed_neighbors = mixed_search.combineNeighborList(mixed_typed_neighbors);
    assert(mixed_neighbors[fluid].size() == 1);
    assert(mixed_neighbors[fluid][0] == wall);
    assert(mixed_typed_neighbors.fluid[fluid].empty());
    assert(mixed_typed_neighbors.wall[fluid].size() == 1);
    assert(mixed_typed_neighbors.wall[fluid][0] == wall);
    assert(mixed_neighbors[wall].size() == 1);
    assert(mixed_neighbors[wall][0] == fluid);
    assert(mixed_typed_neighbors.fluid[wall].size() == 1);
    assert(mixed_typed_neighbors.fluid[wall][0] == fluid);
    assert(mixed_typed_neighbors.wall[wall].empty());
    assert(mixed_neighbors[far_fluid].empty());
    assert(mixed_typed_neighbors.fluid[far_fluid].empty());
    assert(mixed_typed_neighbors.wall[far_fluid].empty());

    mixed_search.updateNeighborCounts(mixed_particles, mixed_typed_neighbors);
    assert(mixed_particles.neighborCounts()[fluid] == 1);
    assert(mixed_particles.fluidNeighborCounts()[fluid] == 0);
    assert(mixed_particles.wallNeighborCounts()[fluid] == 1);
    assert(mixed_particles.neighborCounts()[wall] == 1);
    assert(mixed_particles.fluidNeighborCounts()[wall] == 1);
    assert(mixed_particles.wallNeighborCounts()[wall] == 0);
    assert(mixed_particles.neighborCounts()[far_fluid] == 0);
    assert(mixed_particles.fluidNeighborCounts()[far_fluid] == 0);
    assert(mixed_particles.wallNeighborCounts()[far_fluid] == 0);

    lsmps::ParticleSet large_particles;
    constexpr int fluid_nx = 5;
    constexpr int fluid_ny = 5;
    constexpr int fluid_nz = 5;
    large_particles.reserve(343);
    for (int z = 0; z < fluid_nz; ++z) {
        for (int y = 0; y < fluid_ny; ++y) {
            for (int x = 0; x < fluid_nx; ++x) {
                large_particles.addFluidParticle(
                    {static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)});
            }
        }
    }

    for (int z = -1; z <= fluid_nz; ++z) {
        for (int y = -1; y <= fluid_ny; ++y) {
            for (int x = -1; x <= fluid_nx; ++x) {
                const bool is_fluid_region =
                    x >= 0 && x < fluid_nx && y >= 0 && y < fluid_ny && z >= 0 && z < fluid_nz;
                const bool is_wall_shell =
                    x == -1 || x == fluid_nx || y == -1 || y == fluid_ny || z == -1 || z == fluid_nz;
                if (!is_fluid_region && is_wall_shell) {
                    large_particles.addWallParticle(
                        {static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)});
                }
            }
        }
    }

    assert(large_particles.size() == 343);

    lsmps::NeighborSearch large_search(1.01, {-1.0, -1.0, -1.0});
    const lsmps::TypedNeighborList large_typed_neighbors = large_search.buildTypedNeighborList(large_particles);
    const lsmps::NeighborList large_neighbors = large_search.combineNeighborList(large_typed_neighbors);
    const lsmps::NeighborList large_brute_force_neighbors = buildBruteForceNeighbors(large_particles, 1.01);
    assert(large_neighbors == large_brute_force_neighbors);

    large_search.updateNeighborCounts(large_particles, large_typed_neighbors);
    const std::vector<lsmps::NeighborCountSummary> large_summaries =
        large_search.countNeighborsByType(large_typed_neighbors);

    const std::size_t large_center = gridIndex(2, 2, 2, fluid_nx, fluid_ny);
    assert(large_particles.neighborCounts()[large_center] == 6);
    assert(large_particles.fluidNeighborCounts()[large_center] == 6);
    assert(large_particles.wallNeighborCounts()[large_center] == 0);
    assert(large_typed_neighbors.fluid[large_center].size() == 6);
    assert(large_typed_neighbors.wall[large_center].empty());
    assert(large_summaries[large_center].total == 6);
    assert(large_summaries[large_center].fluid == 6);
    assert(large_summaries[large_center].wall == 0);

    const std::size_t boundary_fluid = gridIndex(0, 2, 2, fluid_nx, fluid_ny);
    assert(large_particles.neighborCounts()[boundary_fluid] == 6);
    assert(large_particles.fluidNeighborCounts()[boundary_fluid] == 5);
    assert(large_particles.wallNeighborCounts()[boundary_fluid] == 1);
    assert(large_typed_neighbors.fluid[boundary_fluid].size() == 5);
    assert(large_typed_neighbors.wall[boundary_fluid].size() == 1);
    assert(large_summaries[boundary_fluid].total == 6);
    assert(large_summaries[boundary_fluid].fluid == 5);
    assert(large_summaries[boundary_fluid].wall == 1);

    std::size_t max_total_neighbors = 0;
    std::size_t particles_with_wall_neighbors = 0;
    for (std::size_t index = 0; index < large_particles.size(); ++index) {
        assert(large_particles.neighborCounts()[index] ==
               large_particles.fluidNeighborCounts()[index] + large_particles.wallNeighborCounts()[index]);
        max_total_neighbors = std::max(max_total_neighbors, large_particles.neighborCounts()[index]);
        if (large_particles.wallNeighborCounts()[index] > 0) {
            ++particles_with_wall_neighbors;
        }
    }
    assert(max_total_neighbors == 6);
    assert(particles_with_wall_neighbors > 0);

    lsmps::CellLinkedList cells({0.0, 0.0, 0.0}, 1.0);
    cells.build(mixed_particles);
    const lsmps::CellIndex origin_cell{0, 0, 0};
    const lsmps::CellIndex positive_x_cell{1, 0, 0};
    const lsmps::CellIndex negative_x_cell{-1, 0, 0};
    assert(cells.cellIndex({0.0, 0.0, 0.0}) == origin_cell);
    assert(cells.cellIndex({1.0, 0.0, 0.0}) == positive_x_cell);
    assert(cells.cellIndex({-0.1, 0.0, 0.0}) == negative_x_cell);
    assert(cells.particlesInCell({0, 0, 0}).size() == 2);

    bool rejected_bad_radius = false;
    try {
        lsmps::NeighborSearch bad_search(0.0);
    } catch (const std::runtime_error&) {
        rejected_bad_radius = true;
    }
    assert(rejected_bad_radius);

    return 0;
}
