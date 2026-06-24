#include "time_integration/time_stepper.hpp"

#include "correction/correction.hpp"
#include "diagnostics/particle_distribution_diagnostics.hpp"
#include "free_surface/free_surface_detector.hpp"
#include "io/vtk_writer.hpp"
#include "lsmps/lsmps_matrices.hpp"
#include "neighbor/neighbor_search.hpp"
#include "particle_shifting/particle_shifter.hpp"
#include "pressure_poisson/pressure_poisson.hpp"
#include "provisional/provisional.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lsmps {
namespace {

using Clock = std::chrono::steady_clock;
using TimePoint = Clock::time_point;

double elapsedSeconds(const TimePoint& begin, const TimePoint& end) {
    return std::chrono::duration<double>(end - begin).count();
}

struct StepTiming {
    double neighbor = 0.0;
    double time_control = 0.0;
    double free_surface = 0.0;
    double lsmps = 0.0;
    double provisional = 0.0;
    double pressure_poisson = 0.0;
    double correction = 0.0;
    double particle_shifting = 0.0;
    double diagnostics = 0.0;
    double vtk_output = 0.0;
    double total = 0.0;
};

void printTiming(std::size_t step, double time, const StepTiming& timing) {
    std::cout << "[timing] step=" << step
              << " time=" << time
              << " neighbor=" << timing.neighbor
              << " time_control=" << timing.time_control
              << " free_surface=" << timing.free_surface
              << " lsmps=" << timing.lsmps
              << " provisional=" << timing.provisional
              << " ppe=" << timing.pressure_poisson
              << " correction=" << timing.correction
              << " particle_shifting=" << timing.particle_shifting
              << " diagnostics=" << timing.diagnostics
              << " vtk=" << timing.vtk_output
              << " total=" << timing.total << '\n';
}

void writeSplitParticleOutput(
    const FileManager& files,
    const ParticleSet& particles,
    std::size_t output_index,
    bool write_wall,
    const std::vector<Vector3>& pressure_gradient = {},
    const ParticleDistributionDiagnostics* distribution = nullptr,
    const std::vector<double>& correction_status = {},
    const ParticleShiftResult* particle_shift = nullptr) {
    const VtkWriter writer;
    VtkBuiltInFieldOptions fluid_fields;
    fluid_fields.particle_type = false;
    fluid_fields.fluid_neighbor_count = false;
    fluid_fields.wall_neighbor_count = false;
    VtkBuiltInFieldOptions wall_fields;
    wall_fields.particle_type = false;
    wall_fields.fluid_state = false;
    wall_fields.fluid_neighbor_count = false;
    wall_fields.wall_neighbor_count = false;
    std::vector<VtkVectorField> vector_fields;
    if (!pressure_gradient.empty()) {
        vector_fields.push_back({"pressure_gradient", pressure_gradient});
    }
    if (particle_shift != nullptr && particle_shift->applied) {
        vector_fields.push_back({"particle_shift", particle_shift->displacement});
    }
    std::vector<VtkScalarField> scalar_fields;
    if (distribution != nullptr) {
        scalar_fields = {
            {"nearest_fluid_distance", distribution->nearest_fluid_distance},
            {"max_fluid_neighbor_distance", distribution->max_fluid_neighbor_distance},
            {"mean_fluid_neighbor_distance", distribution->mean_fluid_neighbor_distance},
            {"number_density", distribution->number_density},
            {"number_density_ratio", distribution->number_density_ratio},
            {"x_positive_neighbor_count", distribution->x_positive_neighbor_count},
            {"x_negative_neighbor_count", distribution->x_negative_neighbor_count},
            {"y_positive_neighbor_count", distribution->y_positive_neighbor_count},
            {"y_negative_neighbor_count", distribution->y_negative_neighbor_count},
            {"z_positive_neighbor_count", distribution->z_positive_neighbor_count},
            {"z_negative_neighbor_count", distribution->z_negative_neighbor_count},
            {"directional_coverage_score", distribution->directional_coverage_score},
            {"geometry_min_eigenvalue", distribution->geometry_min_eigenvalue},
            {"geometry_max_eigenvalue", distribution->geometry_max_eigenvalue},
            {"geometry_condition_number", distribution->geometry_condition_number},
        };
    }
    if (!correction_status.empty()) {
        scalar_fields.push_back({"correction_status", correction_status});
    }
    if (particle_shift != nullptr && particle_shift->applied) {
        scalar_fields.push_back({"particle_shift_magnitude", particle_shift->magnitude});
        scalar_fields.push_back({"particle_shift_limited", particle_shift->limited});
        scalar_fields.push_back({"particle_shift_repulsion_active", particle_shift->repulsion_active});
    }
    writer.writeParticlesByType(
        files.fluidOutputPath(output_index),
        particles,
        ParticleType::Fluid,
        fluid_fields,
        scalar_fields,
        vector_fields);
    if (write_wall) {
        writer.writeParticlesByType(
            files.wallOutputPath(output_index),
            particles,
            ParticleType::Wall,
            wall_fields,
            {},
            {});
    }
}

double maxFluidMagnitude(const ParticleSet& particles, const std::vector<Vector3>& values) {
    double maximum = 0.0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            maximum = std::max(maximum, norm(values[i]));
        }
    }
    return maximum;
}

std::vector<double> correctionStatusField(const CorrectionResult& correction) {
    std::vector<double> result(correction.status.size(), 0.0);
    for (std::size_t i = 0; i < correction.status.size(); ++i) {
        result[i] = static_cast<double>(static_cast<int>(correction.status[i]));
    }
    return result;
}

TimeStepDiagnostics buildDiagnostics(
    const ParticleSet& particles,
    const SimulationConfig& config,
    const SimulationState& state,
    const PressurePoissonResult& pressure,
    const CorrectionResult& correction,
    const TimeStepDecision& time_step_decision,
    double next_dt) {
    TimeStepDiagnostics diagnostics;
    diagnostics.step = state.current_step;
    diagnostics.time = state.current_time;
    diagnostics.dt = time_step_decision.dt;
    diagnostics.next_dt = next_dt;
    diagnostics.cfl_limited_dt = time_step_decision.cfl_limited_dt;
    diagnostics.max_relative_velocity = time_step_decision.max_relative_velocity;
    diagnostics.max_velocity = maxFluidMagnitude(particles, particles.velocities());
    diagnostics.max_displacement = maxFluidMagnitude(particles, correction.displacement);
    diagnostics.cfl_number =
        diagnostics.max_relative_velocity * diagnostics.dt / config.geometry.particle_spacing;
    diagnostics.pressure_converged = pressure.solve_info.converged;
    diagnostics.pressure_iterations = pressure.solve_info.iterations;
    diagnostics.pressure_residual = pressure.solve_info.final_residual_norm;

    bool has_fluid = false;
    std::size_t total_neighbors = 0;
    diagnostics.min_neighbor_count = std::numeric_limits<std::size_t>::max();
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!particles.isFluid(i)) {
            continue;
        }
        has_fluid = true;
        const std::size_t count = particles.neighborCounts()[i];
        diagnostics.min_neighbor_count = std::min(diagnostics.min_neighbor_count, count);
        diagnostics.max_neighbor_count = std::max(diagnostics.max_neighbor_count, count);
        total_neighbors += count;

        switch (particles.fluidStates()[i]) {
        case FluidParticleState::Internal:
            ++diagnostics.internal_count;
            break;
        case FluidParticleState::NearFreeSurface:
            ++diagnostics.near_free_surface_count;
            break;
        case FluidParticleState::FreeSurface:
            ++diagnostics.free_surface_count;
            break;
        case FluidParticleState::Splash:
            ++diagnostics.splash_count;
            break;
        }
    }

    if (has_fluid) {
        const std::size_t fluid_count =
            diagnostics.internal_count + diagnostics.near_free_surface_count +
            diagnostics.free_surface_count + diagnostics.splash_count;
        diagnostics.mean_neighbor_count = static_cast<double>(total_neighbors) / static_cast<double>(fluid_count);
    } else {
        diagnostics.min_neighbor_count = 0;
    }

    return diagnostics;
}

void updateSimulationState(SimulationState& state, const TimeStepDiagnostics& diagnostics) {
    state.max_velocity = diagnostics.max_velocity;
    state.cfl_number = diagnostics.cfl_number;
    state.pressure_residual = diagnostics.pressure_residual;
    state.min_neighbor_count = diagnostics.min_neighbor_count;
    state.max_neighbor_count = diagnostics.max_neighbor_count;
    state.mean_neighbor_count = diagnostics.mean_neighbor_count;
}

}  // namespace

TimeStepper::TimeStepper(SimulationConfig config, TimeStepperOptions options)
    : config_(std::move(config)),
      options_(std::move(options)),
      files_(config_.file),
      time_control_(config_.time, config_.geometry.particle_spacing) {
    const std::vector<std::string> errors = config_.validate();
    if (!errors.empty()) {
        std::string message = "Invalid time-stepper configuration:";
        for (const std::string& error : errors) {
            message += "\n  - " + error;
        }
        throw std::runtime_error(message);
    }
}

const SimulationState& TimeStepper::state() const noexcept {
    return state_;
}

const std::vector<TimeStepDiagnostics>& TimeStepper::history() const noexcept {
    return history_;
}

void TimeStepper::writeCurrentState(const ParticleSet& particles, const std::string& tag) const {
    (void)tag;
    writeSplitParticleOutput(files_, particles, output_index_, true);
}

TimeStepDiagnostics TimeStepper::advanceOneStep(ParticleSet& particles) {
    if (particles.empty()) {
        throw std::runtime_error("TimeStepper requires a non-empty particle set");
    }

    StepTiming timing;
    const TimePoint step_begin = Clock::now();
    TimePoint section_begin = step_begin;

    NeighborSearch neighbor_search(config_.geometry.support_radius, config_.geometry.domain_min);
    TypedNeighborList neighbors = neighbor_search.buildTypedNeighborList(particles);
    neighbor_search.updateNeighborCounts(particles, neighbors);
    const ParticleDistributionDiagnostics distribution_diagnostics =
        computeParticleDistributionDiagnostics(particles, neighbors, config_.geometry, config_.lsmps.kernel_type);
    TimePoint section_end = Clock::now();
    timing.neighbor = elapsedSeconds(section_begin, section_end);

    section_begin = section_end;
    const double max_relative_velocity = computeMaxRelativeVelocity(particles, neighbors);
    const TimeStepDecision time_step_decision =
        time_control_.chooseDt(state_.current_time, max_relative_velocity);
    SimulationConfig step_config = config_;
    step_config.time.dt = time_step_decision.dt;
    section_end = Clock::now();
    timing.time_control = elapsedSeconds(section_begin, section_end);

    section_begin = section_end;
    const FreeSurfaceDetector detector(config_.free_surface, config_.geometry.particle_spacing);
    const FreeSurfaceDiagnostics free_surface = detector.detect(particles, neighbors);
    section_end = Clock::now();
    timing.free_surface = elapsedSeconds(section_begin, section_end);

    section_begin = section_end;
    const LsmpsMatrixSet matrices =
        buildLsmpsMatrices(particles, neighbors, step_config.geometry.support_radius, step_config.lsmps);
    section_end = Clock::now();
    timing.lsmps = elapsedSeconds(section_begin, section_end);

    section_begin = section_end;
    const ProvisionalVelocityCalculator provisional_calculator;
    const ProvisionalVelocityResult provisional =
        provisional_calculator.compute(particles, neighbors, matrices, step_config);
    section_end = Clock::now();
    timing.provisional = elapsedSeconds(section_begin, section_end);

    section_begin = section_end;
    const PressurePoissonAssembler pressure_poisson;
    const PressurePoissonResult pressure =
        pressure_poisson.solve(particles, neighbors, matrices, provisional, step_config);
    if (!pressure.solved || (options_.fail_on_pressure_nonconvergence && !pressure.solve_info.converged)) {
        throw std::runtime_error("PPE solve failed or did not converge during time stepping");
    }
    section_end = Clock::now();
    timing.pressure_poisson = elapsedSeconds(section_begin, section_end);

    section_begin = section_end;
    const PressureCorrectionApplier correction_applier;
    const CorrectionResult correction =
        correction_applier.apply(particles, neighbors, matrices, provisional, pressure, step_config);
    section_end = Clock::now();
    timing.correction = elapsedSeconds(section_begin, section_end);

    section_begin = section_end;
    particles.pressures() = pressure.pressure;
    particles.velocities() = correction.next_velocity;
    particles.positions() = correction.next_position;
    const ParticleShifter particle_shifter;
    const ParticleShiftResult particle_shift =
        particle_shifter.compute(particles, neighbors, config_);
    particle_shifter.apply(particles, particle_shift);
    section_end = Clock::now();
    timing.particle_shifting = elapsedSeconds(section_begin, section_end);

    section_begin = section_end;
    state_.current_step += 1;
    state_.current_time += time_step_decision.dt;
    time_control_.updateAfterStep(time_step_decision);
    const TimeStepDecision next_decision =
        time_control_.chooseDt(state_.current_time, computeMaxRelativeVelocity(particles, neighbors));

    const TimeStepDiagnostics diagnostics =
        buildDiagnostics(particles, config_, state_, pressure, correction, time_step_decision, next_decision.dt);
    updateSimulationState(state_, diagnostics);
    history_.push_back(diagnostics);
    section_end = Clock::now();
    timing.diagnostics = elapsedSeconds(section_begin, section_end);

    if (config_.file.write_outputs && time_control_.shouldWriteOutput(state_.current_time)) {
        section_begin = Clock::now();
        writeSplitParticleOutput(
            files_,
            particles,
            output_index_,
            config_.file.write_wall_each_output,
            correction.pressure_gradient,
            &distribution_diagnostics,
            correctionStatusField(correction),
            &particle_shift);
        section_end = Clock::now();
        timing.vtk_output = elapsedSeconds(section_begin, section_end);
        ++output_index_;
        time_control_.advanceOutputTime(state_.current_time);
        timing.total = elapsedSeconds(step_begin, Clock::now());
        printTiming(state_.current_step, state_.current_time, timing);
    }

    return diagnostics;
}

std::vector<TimeStepDiagnostics> TimeStepper::run(ParticleSet& particles) {
    state_.current_time = config_.time.start_time;
    if (config_.file.write_outputs && config_.file.write_initial_state) {
        writeSplitParticleOutput(files_, particles, output_index_, true);
        ++output_index_;
    }

    while (!time_control_.reachedEndTime(state_.current_time)) {
        advanceOneStep(particles);
    }

    return history_;
}

double TimeStepper::computeMaxRelativeVelocity(
    const ParticleSet& particles,
    const TypedNeighborList& neighbors) const {
    double maximum = 0.0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!particles.isFluid(i)) {
            continue;
        }
        for (const std::size_t j : neighbors.fluid[i]) {
            maximum = std::max(maximum, norm(particles.velocities()[i] - particles.velocities()[j]));
        }
        for (const std::size_t j : neighbors.wall[i]) {
            maximum = std::max(maximum, norm(particles.velocities()[i] - particles.velocities()[j]));
        }
    }
    return maximum;
}

}  // namespace lsmps
