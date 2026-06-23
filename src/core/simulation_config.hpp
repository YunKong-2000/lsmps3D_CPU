#pragma once

#include "core/vector3.hpp"

#include <string>
#include <vector>

namespace lsmps {

struct TimeConfig {
    double dt = 0.001;
    double start_time = 0.0;
    double end_time = 1.0;
    double initial_dt = 0.001;
    double min_dt = 1.0e-6;
    double max_dt = 0.001;
    double cfl_number = 0.2;
    double growth_factor = 1.05;
    double output_interval = 0.01;
};

struct FileConfig {
    std::string input_directory = "cases";
    std::string input_file = "";
    std::string fluid_particle_file = "";
    std::string wall_particle_file = "";
    std::string output_directory = "output/time_stepper";
    std::string output_prefix = "step";
    bool write_initial_state = true;
    bool write_outputs = true;
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
    double screen_radius_factor = 3.1;
    double wall_patch_radius_factor = 0.85;
    double particle_radius_factor = 0.5;
    double open_threshold = 0.18;
    double cone_angle_degrees = 45.0;
    double cone_threshold = 0.62;
    double min_cone_accessible_ratio = 0.40;
    std::size_t cubed_sphere_q = 8;
    std::size_t splash_max_fluid_neighbors = 4;
    double splash_open_threshold = 0.75;
    double near_surface_distance_factor = 1.5;
};

enum class LsmpsKernelType {
    Linear,
};

struct LsmpsConfig {
    std::size_t min_neighbors = 9;
    double eigenvalue_tolerance = 1.0e-12;
    double condition_number_warning = 1.0e8;
    double condition_number_failure = 1.0e12;
    LsmpsKernelType kernel_type = LsmpsKernelType::Linear;
};

struct LinearSolverConfig {
    std::size_t max_iterations = 1000;
    double tolerance = 1.0e-10;
};

struct SimulationConfig {
    TimeConfig time;
    FileConfig file;
    GeometryConfig geometry;
    PhysicalConfig physical;
    FreeSurfaceConfig free_surface;
    LsmpsConfig lsmps;
    LinearSolverConfig linear_solver;

    std::vector<std::string> validate() const;
};

}  // namespace lsmps
