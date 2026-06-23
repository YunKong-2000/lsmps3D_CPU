#include "time_integration/time_stepper.hpp"

#include "correction/correction.hpp"
#include "free_surface/free_surface_detector.hpp"
#include "io/vtk_writer.hpp"
#include "lsmps/lsmps_matrices.hpp"
#include "neighbor/neighbor_search.hpp"
#include "pressure_poisson/pressure_poisson.hpp"
#include "provisional/provisional.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace lsmps {
namespace {

template <typename Enum>
std::vector<double> enumFieldAsDouble(const std::vector<Enum>& values) {
    std::vector<double> result(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        result[i] = static_cast<double>(static_cast<int>(values[i]));
    }
    return result;
}

std::vector<double> intFieldAsDouble(const std::vector<int>& values) {
    std::vector<double> result(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        result[i] = static_cast<double>(values[i]);
    }
    return result;
}

std::vector<double> matrixStatusField(const LsmpsMatrixSet& matrices, bool pressure_neumann) {
    std::vector<double> result(matrices.particles.size(), 0.0);
    for (std::size_t i = 0; i < matrices.particles.size(); ++i) {
        const LsmpsInverseMatrix& matrix =
            pressure_neumann ? matrices.particles[i].pressure_neumann : matrices.particles[i].regular;
        result[i] = static_cast<double>(static_cast<int>(matrix.status));
    }
    return result;
}

std::vector<double> vectorMagnitude(const std::vector<Vector3>& values) {
    std::vector<double> result(values.size(), 0.0);
    for (std::size_t i = 0; i < values.size(); ++i) {
        result[i] = norm(values[i]);
    }
    return result;
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

TimeStepDiagnostics buildDiagnostics(
    const ParticleSet& particles,
    const SimulationConfig& config,
    const SimulationState& state,
    const PressurePoissonResult& pressure,
    const CorrectionResult& correction) {
    TimeStepDiagnostics diagnostics;
    diagnostics.step = state.current_step;
    diagnostics.time = state.current_time;
    diagnostics.dt = config.time.dt;
    diagnostics.max_velocity = maxFluidMagnitude(particles, particles.velocities());
    diagnostics.max_displacement = maxFluidMagnitude(particles, correction.displacement);
    diagnostics.cfl_number = diagnostics.max_velocity * config.time.dt / config.geometry.particle_spacing;
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

std::string formatStepTag(std::size_t output_index, std::size_t step) {
    std::ostringstream stream;
    stream << std::setfill('0') << std::setw(5) << output_index << "_step_" << step;
    return stream.str();
}

}  // namespace

TimeStepper::TimeStepper(SimulationConfig config, TimeStepperOptions options)
    : config_(std::move(config)),
      options_(std::move(options)),
      next_output_time_(config_.time.output_interval) {
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
    const VtkWriter writer;
    writer.writeParticles(outputPath(tag), particles);
}

TimeStepDiagnostics TimeStepper::advanceOneStep(ParticleSet& particles) {
    if (particles.empty()) {
        throw std::runtime_error("TimeStepper requires a non-empty particle set");
    }

    NeighborSearch neighbor_search(config_.geometry.support_radius, config_.geometry.domain_min);
    TypedNeighborList neighbors = neighbor_search.buildTypedNeighborList(particles);
    neighbor_search.updateNeighborCounts(particles, neighbors);

    const FreeSurfaceDetector detector(config_.free_surface, config_.geometry.particle_spacing);
    const FreeSurfaceDiagnostics free_surface = detector.detect(particles, neighbors);

    const LsmpsMatrixSet matrices =
        buildLsmpsMatrices(particles, neighbors, config_.geometry.support_radius, config_.lsmps);

    const ProvisionalVelocityCalculator provisional_calculator;
    const ProvisionalVelocityResult provisional =
        provisional_calculator.compute(particles, neighbors, matrices, config_);

    const PressurePoissonAssembler pressure_poisson;
    const PressurePoissonResult pressure =
        pressure_poisson.solve(particles, neighbors, matrices, provisional, config_);
    if (!pressure.solved || (options_.fail_on_pressure_nonconvergence && !pressure.solve_info.converged)) {
        throw std::runtime_error("PPE solve failed or did not converge during time stepping");
    }

    const PressureCorrectionApplier correction_applier;
    const CorrectionResult correction =
        correction_applier.apply(particles, neighbors, matrices, provisional, pressure, config_);

    particles.pressures() = pressure.pressure;
    particles.velocities() = correction.next_velocity;
    particles.positions() = correction.next_position;

    state_.current_step += 1;
    state_.current_time += config_.time.dt;

    const TimeStepDiagnostics diagnostics = buildDiagnostics(particles, config_, state_, pressure, correction);
    updateSimulationState(state_, diagnostics);
    history_.push_back(diagnostics);

    if (options_.write_outputs && shouldWriteOutput()) {
        const std::string tag = formatStepTag(output_index_, state_.current_step);
        const VtkWriter writer;
        writer.writeParticles(
            outputPath(tag),
            particles,
            {
                {"free_surface_open_ratio", free_surface.open_ratio},
                {"free_surface_cone_ratio", free_surface.cone_ratio},
                {"free_surface_accessible_area_ratio", free_surface.accessible_area_ratio},
                {"free_surface_reason_code", intFieldAsDouble(free_surface.reason_code)},
                {"lsmps_regular_status", matrixStatusField(matrices, false)},
                {"lsmps_pressure_neumann_status", matrixStatusField(matrices, true)},
                {"provisional_status", enumFieldAsDouble(provisional.status)},
                {"correction_status", enumFieldAsDouble(correction.status)},
                {"ppe_rhs", pressure.diagnostics.rhs},
                {"ppe_divergence", pressure.diagnostics.divergence},
                {"ppe_wall_neumann_source", pressure.diagnostics.wall_neumann_source},
                {"pressure_dof", intFieldAsDouble(pressure.diagnostics.pressure_dof)},
                {"pressure_gradient_magnitude", vectorMagnitude(correction.pressure_gradient)},
                {"velocity_correction_magnitude", vectorMagnitude(correction.velocity_correction)},
                {"displacement_magnitude", vectorMagnitude(correction.displacement)},
            },
            {
                {"provisional_velocity", provisional.provisional_velocity},
                {"viscous_acceleration", provisional.viscous_acceleration},
                {"body_acceleration", provisional.body_acceleration},
                {"velocity_delta", provisional.velocity_delta},
                {"pressure_gradient", correction.pressure_gradient},
                {"velocity_correction", correction.velocity_correction},
                {"displacement", correction.displacement},
            });
        ++output_index_;
        while (next_output_time_ <= state_.current_time + 0.5 * config_.time.dt) {
            next_output_time_ += config_.time.output_interval;
        }
    }

    return diagnostics;
}

std::vector<TimeStepDiagnostics> TimeStepper::run(ParticleSet& particles) {
    if (options_.write_outputs && options_.write_initial_state) {
        writeCurrentState(particles, "initial");
    }

    while (state_.current_time + 0.5 * config_.time.dt < config_.time.end_time) {
        advanceOneStep(particles);
    }

    return history_;
}

bool TimeStepper::shouldWriteOutput() const {
    return state_.current_time + 0.5 * config_.time.dt >= next_output_time_ ||
           state_.current_time + 0.5 * config_.time.dt >= config_.time.end_time;
}

std::string TimeStepper::outputPath(const std::string& tag) const {
    const std::filesystem::path directory(options_.output_directory);
    const std::filesystem::path file = options_.output_prefix + "_" + tag + ".vtk";
    return (directory / file).string();
}

}  // namespace lsmps
