#include "io/file_manager.hpp"

#include <cassert>
#include <string>

int main() {
    lsmps::FileConfig config;
    config.input_directory = "cases/demo";
    config.input_file = "initial.vtk";
    config.fluid_particle_file = "fluid.dat";
    config.wall_particle_file = "wall.dat";
    config.output_directory = "output/demo";
    config.output_prefix = "case";

    const lsmps::FileManager manager(config);
    assert(manager.inputPath() == "cases/demo/initial.vtk");
    assert(manager.fluidParticlePath() == "cases/demo/fluid.dat");
    assert(manager.wallParticlePath() == "cases/demo/wall.dat");
    assert(manager.initialOutputPath() == "output/demo/case_initial.vtk");
    assert(manager.stepOutputPath(3, 42) == "output/demo/case_00003.vtk");
    assert(manager.fluidOutputPath(3) == "output/demo/case_fluid_00003.vtk");
    assert(manager.wallOutputPath(3) == "output/demo/case_wall_00003.vtk");

    config.input_file = "/tmp/absolute.vtk";
    const lsmps::FileManager absolute_manager(config);
    assert(absolute_manager.inputPath() == "/tmp/absolute.vtk");

    config.input_file.clear();
    const lsmps::FileManager empty_manager(config);
    assert(empty_manager.inputPath().empty());

    return 0;
}
