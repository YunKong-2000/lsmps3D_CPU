#include "core/particle_set.hpp"
#include "core/particle_types.hpp"
#include "io/vtk_writer.hpp"
#include "neighbor/neighbor_search.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>

namespace {

constexpr int fluid_nx = 50;
constexpr int fluid_ny = 50;
constexpr int fluid_nz = 50;
constexpr double spacing = 1.0;
constexpr double support_radius = 3.1 * spacing;
constexpr double support_radius_squared = support_radius * support_radius;

bool isFluidCoordinate(int x, int y, int z) {
    return x >= 0 && x < fluid_nx && y >= 0 && y < fluid_ny && z >= 0 && z < fluid_nz;
}

bool isWallCoordinate(int x, int y, int z) {
    const bool in_outer_box =
        x >= -1 && x <= fluid_nx && y >= -1 && y <= fluid_ny && z >= -1 && z <= fluid_nz;
    const bool on_wall_shell =
        x == -1 || x == fluid_nx || y == -1 || y == fluid_ny || z == -1 || z == fluid_nz;
    return in_outer_box && on_wall_shell && !isFluidCoordinate(x, y, z);
}

std::size_t fluidGridIndex(int x, int y, int z) {
    return static_cast<std::size_t>((z * fluid_ny + y) * fluid_nx + x);
}

struct ExpectedCounts {
    std::size_t fluid = 0;
    std::size_t wall = 0;
};

ExpectedCounts expectedCountsAt(int x, int y, int z) {
    ExpectedCounts counts;
    const int search_radius = static_cast<int>(std::ceil(support_radius / spacing));

    for (int dz = -search_radius; dz <= search_radius; ++dz) {
        for (int dy = -search_radius; dy <= search_radius; ++dy) {
            for (int dx = -search_radius; dx <= search_radius; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) {
                    continue;
                }

                const double distance_squared = spacing * spacing * static_cast<double>(dx * dx + dy * dy + dz * dz);
                if (distance_squared > support_radius_squared) {
                    continue;
                }

                const int neighbor_x = x + dx;
                const int neighbor_y = y + dy;
                const int neighbor_z = z + dz;
                if (isFluidCoordinate(neighbor_x, neighbor_y, neighbor_z)) {
                    ++counts.fluid;
                } else if (isWallCoordinate(neighbor_x, neighbor_y, neighbor_z)) {
                    ++counts.wall;
                }
            }
        }
    }

    return counts;
}

}  // namespace

int main() {
    constexpr std::size_t fluid_particle_count =
        static_cast<std::size_t>(fluid_nx) * fluid_ny * fluid_nz;
    constexpr std::size_t wall_particle_count =
        static_cast<std::size_t>(fluid_nx + 2) * (fluid_ny + 2) * (fluid_nz + 2) - fluid_particle_count;
    constexpr std::size_t particle_count = fluid_particle_count + wall_particle_count;

    lsmps::ParticleSet particles;
    particles.reserve(particle_count);

    for (int z = 0; z < fluid_nz; ++z) {
        for (int y = 0; y < fluid_ny; ++y) {
            for (int x = 0; x < fluid_nx; ++x) {
                particles.addFluidParticle(
                    {spacing * x, spacing * y, spacing * z},
                    {},
                    lsmps::FluidParticleState::Internal,
                    0.0,
                    1000.0);
            }
        }
    }

    for (int z = -1; z <= fluid_nz; ++z) {
        for (int y = -1; y <= fluid_ny; ++y) {
            for (int x = -1; x <= fluid_nx; ++x) {
                if (isWallCoordinate(x, y, z)) {
                    particles.addWallParticle({spacing * x, spacing * y, spacing * z});
                }
            }
        }
    }

    assert(particles.size() == particle_count);
    assert(fluid_particle_count == 125000);
    assert(wall_particle_count == 15608);

    lsmps::NeighborSearch search(support_radius, {-spacing, -spacing, -spacing});
    const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
    search.updateNeighborCounts(particles, neighbors);

    std::size_t min_fluid_neighbor_count = particles.size();
    std::size_t max_fluid_neighbor_count = 0;
    std::size_t min_wall_neighbor_count = particles.size();
    std::size_t max_wall_neighbor_count = 0;
    std::size_t fluid_particles_with_wall_neighbors = 0;
    std::size_t total_fluid_neighbor_links = 0;
    std::size_t total_wall_neighbor_links = 0;

    for (int z = 0; z < fluid_nz; ++z) {
        for (int y = 0; y < fluid_ny; ++y) {
            for (int x = 0; x < fluid_nx; ++x) {
                const std::size_t index = fluidGridIndex(x, y, z);
                const ExpectedCounts expected = expectedCountsAt(x, y, z);
                const std::size_t expected_total = expected.fluid + expected.wall;

                assert(particles.types()[index] == lsmps::ParticleType::Fluid);
                assert(neighbors.fluid[index].size() == expected.fluid);
                assert(neighbors.wall[index].size() == expected.wall);
                assert(neighbors.fluid[index].size() + neighbors.wall[index].size() == expected_total);
                assert(particles.neighborCounts()[index] == expected_total);
                assert(particles.fluidNeighborCounts()[index] == expected.fluid);
                assert(particles.wallNeighborCounts()[index] == expected.wall);

                min_fluid_neighbor_count = std::min(min_fluid_neighbor_count, expected.fluid);
                max_fluid_neighbor_count = std::max(max_fluid_neighbor_count, expected.fluid);
                min_wall_neighbor_count = std::min(min_wall_neighbor_count, expected.wall);
                max_wall_neighbor_count = std::max(max_wall_neighbor_count, expected.wall);
                total_fluid_neighbor_links += expected.fluid;
                total_wall_neighbor_links += expected.wall;
                if (expected.wall > 0) {
                    ++fluid_particles_with_wall_neighbors;
                }
            }
        }
    }

    const ExpectedCounts center_expected = expectedCountsAt(25, 25, 25);
    assert(center_expected.fluid == 122);
    assert(center_expected.wall == 0);
    assert(neighbors.fluid[fluidGridIndex(25, 25, 25)].size() == center_expected.fluid);
    assert(neighbors.wall[fluidGridIndex(25, 25, 25)].empty());
    assert(particles.fluidNeighborCounts()[fluidGridIndex(25, 25, 25)] == center_expected.fluid);
    assert(particles.wallNeighborCounts()[fluidGridIndex(25, 25, 25)] == center_expected.wall);

    const ExpectedCounts corner_expected = expectedCountsAt(0, 0, 0);
    assert(corner_expected.wall > 0);
    assert(neighbors.fluid[fluidGridIndex(0, 0, 0)].size() == corner_expected.fluid);
    assert(neighbors.wall[fluidGridIndex(0, 0, 0)].size() == corner_expected.wall);
    assert(particles.fluidNeighborCounts()[fluidGridIndex(0, 0, 0)] == corner_expected.fluid);
    assert(particles.wallNeighborCounts()[fluidGridIndex(0, 0, 0)] == corner_expected.wall);

    const lsmps::VtkWriter writer;
    const std::string output_path = "output/neighbor_search_50x50x50.vtk";
    writer.writeParticles(output_path, particles);

    std::cout << "Wrote large neighbor-search VTK output: " << output_path << '\n';
    std::cout << "Fluid particle count: " << fluid_particle_count << '\n';
    std::cout << "Wall particle count: " << wall_particle_count << '\n';
    std::cout << "Support radius: " << support_radius << '\n';
    std::cout << "Fluid neighbor count range on fluid particles: " << min_fluid_neighbor_count << " to "
              << max_fluid_neighbor_count << '\n';
    std::cout << "Wall neighbor count range on fluid particles: " << min_wall_neighbor_count << " to "
              << max_wall_neighbor_count << '\n';
    std::cout << "Center fluid particle counts: fluid=" << center_expected.fluid
              << ", wall=" << center_expected.wall << '\n';
    std::cout << "Corner fluid particle counts: fluid=" << corner_expected.fluid
              << ", wall=" << corner_expected.wall << '\n';
    std::cout << "Fluid particles with wall neighbors: " << fluid_particles_with_wall_neighbors << '\n';
    std::cout << "Total fluid-neighbor links on fluid particles: " << total_fluid_neighbor_links << '\n';
    std::cout << "Total wall-neighbor links on fluid particles: " << total_wall_neighbor_links << '\n';

    return 0;
}
