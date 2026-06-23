#include "core/particle_set.hpp"
#include "io/particle_file_io.hpp"

#include <cassert>

int main() {
    lsmps::ParticleSet source;
    source.addFluidParticle({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, lsmps::FluidParticleState::Internal, 2.0, 1000.0);
    source.addWallParticle({0.0, -1.0, 0.0}, {}, 0.0, 1000.0, {0.0, 1.0, 0.0});

    const lsmps::ParticleFileWriter writer;
    writer.writeFluidParticles("output/tests/particle_file_io_fluid.dat", source);
    writer.writeWallParticles("output/tests/particle_file_io_wall.dat", source);

    lsmps::ParticleSet loaded;
    const lsmps::ParticleFileReader reader;
    reader.readFluidParticles("output/tests/particle_file_io_fluid.dat", loaded);
    reader.readWallParticles("output/tests/particle_file_io_wall.dat", loaded);

    assert(loaded.size() == 2);
    assert(loaded.isFluid(0));
    assert(loaded.isWall(1));
    assert(loaded.positions()[0].x == 0.0);
    assert(loaded.velocities()[0].x == 1.0);
    assert(loaded.pressures()[0] == 2.0);
    assert(loaded.densities()[0] == 1000.0);
    assert(loaded.wallNormals()[1].y == 1.0);

    return 0;
}
