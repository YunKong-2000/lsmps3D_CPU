#include "core/particle_set.hpp"
#include "core/vector3.hpp"
#include "io/particle_file_io.hpp"
#include "io/vtk_writer.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <string>

namespace {

constexpr double length = 1.0;
constexpr double height = 1.0;
constexpr double width = 1.0;
constexpr double water_height = 0.85;
constexpr double spacing = 0.05;
constexpr double density = 1000.0;

int pointCount(double extent) {
    return static_cast<int>(std::round(extent / spacing)) + 1;
}

double coordinate(int index) {
    return static_cast<double>(index) * spacing;
}

void addFluid(lsmps::ParticleSet& particles) {
    const int nx = pointCount(length);
    const int ny = pointCount(water_height);
    const int nz = pointCount(width);

    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                particles.addFluidParticle(
                    {coordinate(x), coordinate(y), coordinate(z)},
                    {},
                    lsmps::FluidParticleState::Internal,
                    0.0,
                    density);
            }
        }
    }
}

void addWalls(lsmps::ParticleSet& particles) {
    const int nx = pointCount(length);
    const int ny = pointCount(height);
    const int nz = pointCount(width);

    for (int z = -1; z <= nz; ++z) {
        for (int y = -1; y < ny; ++y) {
            for (int x = -1; x <= nx; ++x) {
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
                if (z == -1) {
                    normal += {0.0, 0.0, 1.0};
                    is_wall = true;
                }
                if (z == nz) {
                    normal += {0.0, 0.0, -1.0};
                    is_wall = true;
                }
                if (!is_wall) {
                    continue;
                }

                const double length_normal = lsmps::norm(normal);
                if (length_normal > 0.0) {
                    normal /= length_normal;
                }
                particles.addWallParticle(
                    {coordinate(x), coordinate(y), coordinate(z)},
                    {},
                    0.0,
                    density,
                    normal);
            }
        }
    }
}

}  // namespace

int main() {
    try {
        lsmps::ParticleSet particles;
        particles.reserve(30000);
        addFluid(particles);
        addWalls(particles);

        const lsmps::ParticleFileWriter writer;
        writer.writeFluidParticles("cases/hydrostatic_box/fluid_particles.dat", particles);
        writer.writeWallParticles("cases/hydrostatic_box/wall_particles.dat", particles);

        const lsmps::VtkWriter vtk_writer;
        vtk_writer.writeParticles("output/hydrostatic_box/preprocess_particles.vtk", particles);

        std::cout << "Generated hydrostatic box particles\n";
        std::cout << "  tank = 1.0m x 1.0m x 1.0m, water height = 0.85m\n";
        std::cout << "  fluid + wall total = " << particles.size() << '\n';
        std::cout << "  fluid file = cases/hydrostatic_box/fluid_particles.dat\n";
        std::cout << "  wall file = cases/hydrostatic_box/wall_particles.dat\n";
        std::cout << "  preview vtk = output/hydrostatic_box/preprocess_particles.vtk\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
