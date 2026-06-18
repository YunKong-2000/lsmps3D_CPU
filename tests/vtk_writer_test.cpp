#include "core/particle_set.hpp"
#include "io/vtk_writer.hpp"

#include <cassert>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

std::string readFile(const std::string& path) {
    std::ifstream file(path);
    assert(file);
    std::ostringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

bool contains(const std::string& text, const std::string& pattern) {
    return text.find(pattern) != std::string::npos;
}

}  // namespace

int main() {
    lsmps::ParticleSet particles;
    particles.addFluidParticle(
        {0.0, 0.0, 0.0},
        {1.0, 0.0, 0.0},
        lsmps::FluidParticleState::FreeSurface,
        10.0,
        1000.0);
    particles.addWallParticle({1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 20.0, 1000.0);
    particles.neighborCounts()[0] = 12;
    particles.neighborCounts()[1] = 8;
    particles.fluidNeighborCounts()[0] = 9;
    particles.fluidNeighborCounts()[1] = 5;
    particles.wallNeighborCounts()[0] = 3;
    particles.wallNeighborCounts()[1] = 3;

    const std::string path = "vtk_writer_test_output.vtk";
    const lsmps::VtkWriter writer;
    writer.writeParticles(
        path,
        particles,
        {{"diagnostic_scalar", {0.5, 1.5}}},
        {{"diagnostic_vector", {{0.0, 0.0, 1.0}, {1.0, 1.0, 0.0}}}});

    const std::string output = readFile(path);
    assert(contains(output, "# vtk DataFile Version 3.0\n"));
    assert(contains(output, "DATASET POLYDATA\n"));
    assert(contains(output, "POINTS 2 double\n"));
    assert(contains(output, "VERTICES 2 4\n"));
    assert(contains(output, "POINT_DATA 2\n"));
    assert(contains(output, "VECTORS velocity double\n"));
    assert(contains(output, "SCALARS pressure double 1\n"));
    assert(contains(output, "SCALARS particle_type int 1\nLOOKUP_TABLE default\n0\n1\n"));
    assert(contains(output, "SCALARS fluid_state int 1\nLOOKUP_TABLE default\n2\n0\n"));
    assert(contains(output, "SCALARS neighbor_count int 1\nLOOKUP_TABLE default\n12\n8\n"));
    assert(contains(output, "SCALARS fluid_neighbor_count int 1\nLOOKUP_TABLE default\n9\n5\n"));
    assert(contains(output, "SCALARS wall_neighbor_count int 1\nLOOKUP_TABLE default\n3\n3\n"));
    assert(contains(output, "SCALARS diagnostic_scalar double 1\n"));
    assert(contains(output, "VECTORS diagnostic_vector double\n"));

    bool rejected_bad_size = false;
    try {
        writer.writeParticles("bad_vtk_writer_test_output.vtk", particles, {{"bad_scalar", {1.0}}});
    } catch (const std::runtime_error&) {
        rejected_bad_size = true;
    }
    assert(rejected_bad_size);

    bool rejected_bad_name = false;
    try {
        writer.writeParticles("bad_vtk_writer_test_output.vtk", particles, {{"bad scalar", {1.0, 2.0}}});
    } catch (const std::runtime_error&) {
        rejected_bad_name = true;
    }
    assert(rejected_bad_name);

    return 0;
}
