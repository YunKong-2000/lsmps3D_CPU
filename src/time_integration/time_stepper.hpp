#pragma once

#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "core/simulation_state.hpp"
#include "io/file_manager.hpp"
#include "neighbor/neighbor_search.hpp"
#include "time_integration/time_step_controller.hpp"

#include <cstddef>
#include <vector>

namespace lsmps {

struct TimeStepperOptions {
    bool fail_on_pressure_nonconvergence = true;
};

struct TimeStepDiagnostics {
    std::size_t step = 0;
    double time = 0.0;
    double dt = 0.0;
    double next_dt = 0.0;
    double cfl_limited_dt = 0.0;
    double max_relative_velocity = 0.0;
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
    FileManager files_;
    TimeStepController time_control_;
    SimulationState state_;
    std::vector<TimeStepDiagnostics> history_;
    std::size_t output_index_ = 0;

    double computeMaxRelativeVelocity(
        const ParticleSet& particles,
        const TypedNeighborList& neighbors) const;
};

}  // namespace lsmps
