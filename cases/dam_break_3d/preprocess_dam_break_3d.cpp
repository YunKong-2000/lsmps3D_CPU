#include "core/particle_set.hpp"
#include "core/vector3.hpp"
#include "io/particle_file_io.hpp"
#include "io/vtk_writer.hpp"

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>

namespace {

constexpr double tank_length = 3.2;
constexpr double tank_width = 0.5;
constexpr double tank_height = 1.0;
constexpr double water_length = 1.2;
constexpr double water_height = 0.6;
constexpr double spacing = 0.05;
constexpr double density = 1000.0;

int pointCount(double extent) {
    return static_cast<int>(std::floor(extent / spacing + 1.0e-12)) + 1;
}

double coordinate(int index) {
    return static_cast<double>(index) * spacing;
}

double wallCoordinate(int index, int count, double extent) {
    if (index == -1) {
        return -spacing;
    }
    if (index == count) {
        return extent + spacing;
    }
    return coordinate(index);
}

void addFluid(lsmps::ParticleSet& particles) {
    const int nx = pointCount(water_length);
    const int ny = pointCount(water_height);
    const int nz = pointCount(tank_width);

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
    const int nx = pointCount(tank_length);
    const int ny = pointCount(tank_height);
    const int nz = pointCount(tank_width);

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

                const double normal_length = lsmps::norm(normal);
                if (normal_length > 0.0) {
                    normal /= normal_length;
                }

                particles.addWallParticle(
                    {
                        wallCoordinate(x, nx, tank_length),
                        wallCoordinate(y, ny, tank_height),
                        wallCoordinate(z, nz, tank_width),
                    },
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
        particles.reserve(1000);
        addFluid(particles);
        addWalls(particles);

        const lsmps::ParticleFileWriter writer;
        writer.writeFluidParticles("cases/dam_break_3d/fluid_particles.dat", particles);
        writer.writeWallParticles("cases/dam_break_3d/wall_particles.dat", particles);

        const lsmps::VtkWriter vtk_writer;
        vtk_writer.writeParticles("output/dam_break_3d/preprocess_particles.vtk", particles);

        std::size_t fluid_count = 0;
        std::size_t wall_count = 0;
        for (std::size_t i = 0; i < particles.size(); ++i) {
            if (particles.isFluid(i)) {
                ++fluid_count;
            } else if (particles.isWall(i)) {
                ++wall_count;
            }
        }

        std::cout << "Generated 3D dam-break particles\n";
        std::cout << "  tank = 3.2m x 0.5m x 1.0m (length x width x height)\n";
        std::cout << "  water block = 1.2m x 0.5m x 0.6m\n";
        std::cout << "  spacing = 0.05m\n";
        std::cout << "  fluid particles = " << fluid_count << '\n';
        std::cout << "  wall particles = " << wall_count << '\n';
        std::cout << "  fluid file = cases/dam_break_3d/fluid_particles.dat\n";
        std::cout << "  wall file = cases/dam_break_3d/wall_particles.dat\n";
        std::cout << "  preview vtk = output/dam_break_3d/preprocess_particles.vtk\n";
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
