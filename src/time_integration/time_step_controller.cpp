#include "time_integration/time_step_controller.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace lsmps {
namespace {

constexpr double time_epsilon = 1.0e-12;

}  // namespace

TimeStepController::TimeStepController(TimeStepControlConfig config, double particle_spacing)
    : config_(std::move(config)),
      particle_spacing_(particle_spacing),
      current_dt_(config_.initial_dt),
      next_output_time_(config_.start_time + config_.output_interval) {
    if (particle_spacing_ <= 0.0) {
        throw std::runtime_error("TimeStepController particle spacing must be positive");
    }
}

const TimeStepControlConfig& TimeStepController::config() const noexcept {
    return config_;
}

double TimeStepController::currentDt() const noexcept {
    return current_dt_;
}

TimeStepDecision TimeStepController::chooseDt(double current_time, double max_relative_velocity) {
    TimeStepDecision decision;
    decision.max_relative_velocity = max_relative_velocity;
    decision.cfl_limited_dt = cflLimit(max_relative_velocity);

    double dt = current_dt_;
    if (dt < decision.cfl_limited_dt && dt < config_.max_dt) {
        dt *= config_.growth_factor;
    }

    if (dt > config_.max_dt) {
        dt = config_.max_dt;
        decision.limited_by_max_dt = true;
    }
    if (dt > decision.cfl_limited_dt) {
        dt = decision.cfl_limited_dt;
        decision.limited_by_cfl = true;
    }
    if (dt < config_.min_dt) {
        dt = config_.min_dt;
        decision.limited_by_min_dt = true;
    }

    const double time_to_output = next_output_time_ - current_time;
    if (time_to_output > time_epsilon && dt > time_to_output) {
        dt = time_to_output;
        decision.limited_by_output_time = true;
    }

    const double time_to_end = config_.end_time - current_time;
    if (time_to_end > time_epsilon && dt > time_to_end) {
        dt = time_to_end;
        decision.limited_by_end_time = true;
    }

    decision.dt = dt;
    return decision;
}

void TimeStepController::updateAfterStep(const TimeStepDecision& decision) {
    current_dt_ = decision.dt;
}

bool TimeStepController::reachedEndTime(double current_time) const {
    return current_time >= config_.end_time - time_epsilon;
}

bool TimeStepController::shouldWriteOutput(double current_time) const {
    return current_time >= next_output_time_ - time_epsilon ||
           current_time >= config_.end_time - time_epsilon;
}

void TimeStepController::advanceOutputTime(double current_time) {
    while (next_output_time_ <= current_time + time_epsilon) {
        next_output_time_ += config_.output_interval;
    }
}

double TimeStepController::cflLimit(double max_relative_velocity) const {
    if (max_relative_velocity <= std::numeric_limits<double>::epsilon()) {
        return config_.max_dt;
    }
    return config_.cfl_number * particle_spacing_ / max_relative_velocity;
}

}  // namespace lsmps
