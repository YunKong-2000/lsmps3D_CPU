#include "core/particle_set.hpp"
#include "core/particle_types.hpp"
#include "core/simulation_config.hpp"
#include "free_surface/free_surface_detector.hpp"
#include "io/vtk_writer.hpp"
#include "neighbor/neighbor_search.hpp"

#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

namespace {

constexpr double spacing = 1.0;

std::size_t gridIndex(int x, int y, int z, int nx, int ny) {
    return static_cast<std::size_t>((z * ny + y) * nx + x);
}

lsmps::FreeSurfaceConfig testConfig() {
    lsmps::FreeSurfaceConfig config;
    config.screen_radius_factor = 3.1;
    config.wall_patch_radius_factor = 0.85;
    config.particle_radius_factor = 0.5;
    config.open_threshold = 0.18;
    config.cone_angle_degrees = 45.0;
    config.cone_threshold = 0.62;
    config.min_cone_accessible_ratio = 0.40;
    config.cubed_sphere_q = 8;
    config.splash_max_fluid_neighbors = 4;
    config.splash_open_threshold = 0.75;
    config.near_surface_distance_factor = 1.5;
    return config;
}

void addFluidBlock(lsmps::ParticleSet& particles, int nx, int ny, int nz) {
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                particles.addFluidParticle({spacing * x, spacing * y, spacing * z});
            }
        }
    }
}

void addBoxWalls(lsmps::ParticleSet& particles, int nx, int ny, int nz, bool include_top) {
    for (int z = -1; z <= nz; ++z) {
        for (int y = -1; y <= ny; ++y) {
            for (int x = -1; x <= nx; ++x) {
                const bool in_fluid = x >= 0 && x < nx && y >= 0 && y < ny && z >= 0 && z < nz;
                if (in_fluid) {
                    continue;
                }

                lsmps::Vector3 normal;
                bool is_wall = false;
                if (x == -1) {
                    normal += {1.0, 0.0, 0.0};
                    is_wall = true;
                }
                if (x == nx) {
                    normal += {-1.0, 0.0, 0.0};
                    is_wall = true;
                }
                if (y == -1) {
                    normal += {0.0, 1.0, 0.0};
                    is_wall = true;
                }
                if (y == ny) {
                    normal += {0.0, -1.0, 0.0};
                    is_wall = true;
                }
                if (z == -1) {
                    normal += {0.0, 0.0, 1.0};
                    is_wall = true;
                }
                if (include_top && z == nz) {
                    normal += {0.0, 0.0, -1.0};
                    is_wall = true;
                }

                if (is_wall) {
                    const double length = lsmps::norm(normal);
                    if (length > 0.0) {
                        normal /= length;
                    }
                    particles.addWallParticle({spacing * x, spacing * y, spacing * z}, {}, 0.0, 0.0, normal);
                }
            }
        }
    }
}

std::vector<double> toDoubleVector(const std::vector<int>& values) {
    std::vector<double> result;
    result.reserve(values.size());
    for (const int value : values) {
        result.push_back(static_cast<double>(value));
    }
    return result;
}

void writeFreeSurfaceVtk(
    const std::string& path,
    const lsmps::ParticleSet& particles,
    const lsmps::FreeSurfaceDiagnostics& diagnostics) {
    const lsmps::VtkWriter writer;
    writer.writeParticles(
        path,
        particles,
        {
            {"free_surface_open_ratio", diagnostics.open_ratio},
            {"free_surface_cone_ratio", diagnostics.cone_ratio},
            {"free_surface_accessible_area_ratio", diagnostics.accessible_area_ratio},
            {"free_surface_reason_code", toDoubleVector(diagnostics.reason_code)},
        });
}

}  // namespace

int main() {
    assert(static_cast<int>(lsmps::FluidParticleState::Internal) == 0);
    assert(static_cast<int>(lsmps::FluidParticleState::NearFreeSurface) == 1);
    assert(static_cast<int>(lsmps::FluidParticleState::FreeSurface) == 2);
    assert(static_cast<int>(lsmps::FluidParticleState::Splash) == 3);

    const lsmps::FreeSurfaceConfig config = testConfig();
    const lsmps::FreeSurfaceDetector detector(config, spacing);

    {
        lsmps::ParticleSet particles;
        constexpr int nx = 7;
        constexpr int ny = 7;
        constexpr int nz = 5;
        particles.reserve(nx * ny * nz);
        addFluidBlock(particles, nx, ny, nz);

        lsmps::NeighborSearch search(config.screen_radius_factor * spacing, {0.0, 0.0, 0.0});
        const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
        search.updateNeighborCounts(particles, neighbors);
        const lsmps::FreeSurfaceDiagnostics diagnostics = detector.detect(particles, neighbors);
        writeFreeSurfaceVtk("output/free_surface_open_top.vtk", particles, diagnostics);

        const std::size_t top_center = gridIndex(3, 3, nz - 1, nx, ny);
        const std::size_t interior = gridIndex(3, 3, 2, nx, ny);
        assert(particles.fluidStates()[top_center] == lsmps::FluidParticleState::FreeSurface);
        assert(diagnostics.open_ratio[top_center] >= config.open_threshold);
        assert(diagnostics.cone_ratio[top_center] >= config.cone_threshold);
        assert(particles.fluidStates()[interior] != lsmps::FluidParticleState::FreeSurface);
    }

    {
        lsmps::ParticleSet particles;
        constexpr int nx = 7;
        constexpr int ny = 7;
        constexpr int nz = 7;
        particles.reserve(729);
        addFluidBlock(particles, nx, ny, nz);
        addBoxWalls(particles, nx, ny, nz, true);

        lsmps::NeighborSearch search(config.screen_radius_factor * spacing, {-1.0, -1.0, -1.0});
        const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
        search.updateNeighborCounts(particles, neighbors);
        const lsmps::FreeSurfaceDiagnostics diagnostics = detector.detect(particles, neighbors);
        writeFreeSurfaceVtk("output/free_surface_closed_box.vtk", particles, diagnostics);

        const std::size_t center = gridIndex(3, 3, 3, nx, ny);
        const std::size_t near_wall = gridIndex(0, 3, 3, nx, ny);
        assert(particles.fluidStates()[center] != lsmps::FluidParticleState::FreeSurface);
        assert(particles.fluidStates()[near_wall] != lsmps::FluidParticleState::FreeSurface);
    }

    {
        lsmps::ParticleSet particles;
        constexpr int nx = 7;
        constexpr int ny = 7;
        constexpr int nz = 5;
        particles.reserve(600);
        addFluidBlock(particles, nx, ny, nz);
        addBoxWalls(particles, nx, ny, nz, false);

        lsmps::NeighborSearch search(config.screen_radius_factor * spacing, {-1.0, -1.0, -1.0});
        const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
        search.updateNeighborCounts(particles, neighbors);
        const lsmps::FreeSurfaceDiagnostics diagnostics = detector.detect(particles, neighbors);
        writeFreeSurfaceVtk("output/free_surface_wall_contact.vtk", particles, diagnostics);

        const std::size_t top_near_wall = gridIndex(0, 3, nz - 1, nx, ny);
        const std::size_t side_internal = gridIndex(0, 3, 2, nx, ny);
        assert(particles.fluidStates()[top_near_wall] == lsmps::FluidParticleState::FreeSurface);
        assert(particles.fluidStates()[side_internal] != lsmps::FluidParticleState::FreeSurface);
    }

    {
        lsmps::ParticleSet particles;
        const std::size_t splash = particles.addFluidParticle({0.0, 0.0, 0.0});
        particles.addFluidParticle({spacing, 0.0, 0.0});
        particles.addFluidParticle({0.0, spacing, 0.0});

        lsmps::NeighborSearch search(config.screen_radius_factor * spacing, {0.0, 0.0, 0.0});
        const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
        search.updateNeighborCounts(particles, neighbors);
        const lsmps::FreeSurfaceDiagnostics diagnostics = detector.detect(particles, neighbors);
        writeFreeSurfaceVtk("output/free_surface_splash.vtk", particles, diagnostics);

        assert(particles.fluidStates()[splash] == lsmps::FluidParticleState::Splash);
        assert(diagnostics.reason_code[splash] == 4);
        assert(diagnostics.open_ratio[splash] >= config.splash_open_threshold);
    }

    {
        lsmps::ParticleSet particles;
        constexpr int nx = 7;
        constexpr int ny = 7;
        constexpr int nz = 5;
        particles.reserve(nx * ny * nz);
        addFluidBlock(particles, nx, ny, nz);

        lsmps::FreeSurfaceConfig near_config = config;
        near_config.near_surface_distance_factor = 0.5;
        const lsmps::FreeSurfaceDetector near_detector(near_config, spacing);
        lsmps::NeighborSearch search(near_config.screen_radius_factor * spacing, {0.0, 0.0, 0.0});
        const lsmps::TypedNeighborList neighbors = search.buildTypedNeighborList(particles);
        near_detector.detect(particles, neighbors);

        const std::size_t top_center = gridIndex(3, 3, nz - 1, nx, ny);
        const std::size_t one_layer_below_top = gridIndex(3, 3, nz - 2, nx, ny);
        assert(particles.fluidStates()[top_center] == lsmps::FluidParticleState::FreeSurface);
        assert(particles.fluidStates()[one_layer_below_top] != lsmps::FluidParticleState::NearFreeSurface);
    }

    return 0;
}
