#pragma once

#include "core/vector3.hpp"

#include <string>
#include <vector>

namespace lsmps {

struct TimeConfig {
    double dt = 0.001;
    double end_time = 1.0;
    double output_interval = 0.01;
};

struct GeometryConfig {
    double particle_spacing = 0.02;
    double support_radius = 0.06;
    Vector3 domain_min = {0.0, 0.0, 0.0};
    Vector3 domain_max = {1.0, 1.0, 1.0};
    Vector3 fluid_min = {0.0, 0.0, 0.0};
    Vector3 fluid_max = {0.5, 1.0, 0.5};
};

struct PhysicalConfig {
    double density = 1000.0;
    double viscosity = 1.0e-6;
    Vector3 gravity = {0.0, -9.81, 0.0};
};

struct FreeSurfaceConfig {
    double neighbor_count_ratio = 0.85;
    double number_density_ratio = 0.90;
    std::size_t near_surface_layers = 1;
};

struct LinearSolverConfig {
    std::size_t max_iterations = 1000;
    double tolerance = 1.0e-10;
};

struct SimulationConfig {
    TimeConfig time;
    GeometryConfig geometry;
    PhysicalConfig physical;
    FreeSurfaceConfig free_surface;
    LinearSolverConfig linear_solver;

    std::vector<std::string> validate() const;
};

}  // namespace lsmps
