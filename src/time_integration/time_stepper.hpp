#pragma once

#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "core/simulation_state.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace lsmps {

struct TimeStepperOptions {
    std::string output_directory = "output/time_stepper";
    std::string output_prefix = "step";
    bool write_initial_state = true;
    bool write_outputs = true;
    bool fail_on_pressure_nonconvergence = true;
};

struct TimeStepDiagnostics {
    std::size_t step = 0;
    double time = 0.0;
    double dt = 0.0;
    double max_velocity = 0.0;
    double max_displacement = 0.0;
    double cfl_number = 0.0;
    std::size_t min_neighbor_count = 0;
    std::size_t max_neighbor_count = 0;
    double mean_neighbor_count = 0.0;
    std::size_t internal_count = 0;
    std::size_t near_free_surface_count = 0;
    std::size_t free_surface_count = 0;
    std::size_t splash_count = 0;
    bool pressure_converged = false;
    int pressure_iterations = 0;
    double pressure_residual = 0.0;
};

class TimeStepper {
public:
    explicit TimeStepper(SimulationConfig config, TimeStepperOptions options = {});

    const SimulationState& state() const noexcept;
    const std::vector<TimeStepDiagnostics>& history() const noexcept;

    void writeCurrentState(const ParticleSet& particles, const std::string& tag) const;
    TimeStepDiagnostics advanceOneStep(ParticleSet& particles);
    std::vector<TimeStepDiagnostics> run(ParticleSet& particles);

private:
    SimulationConfig config_;
    TimeStepperOptions options_;
    SimulationState state_;
    std::vector<TimeStepDiagnostics> history_;
    double next_output_time_ = 0.0;
    std::size_t output_index_ = 0;

    bool shouldWriteOutput() const;
    std::string outputPath(const std::string& tag) const;
};

}  // namespace lsmps
