#pragma once

#include "core/simulation_config.hpp"

namespace lsmps {

struct TimeStepDecision {
    double dt = 0.0;
    double cfl_limited_dt = 0.0;
    double max_relative_velocity = 0.0;
    bool limited_by_cfl = false;
    bool limited_by_max_dt = false;
    bool limited_by_min_dt = false;
    bool limited_by_output_time = false;
    bool limited_by_end_time = false;
};

class TimeStepController {
public:
    explicit TimeStepController(TimeConfig config = {}, double particle_spacing = 1.0);

    const TimeConfig& config() const noexcept;
    double currentDt() const noexcept;

    TimeStepDecision chooseDt(double current_time, double max_relative_velocity);
    void updateAfterStep(const TimeStepDecision& decision);
    bool reachedEndTime(double current_time) const;
    bool shouldWriteOutput(double current_time) const;
    void advanceOutputTime(double current_time);

private:
    TimeConfig config_;
    double particle_spacing_ = 1.0;
    double current_dt_ = 0.0;
    double next_output_time_ = 0.0;

    double cflLimit(double max_relative_velocity) const;
};

}  // namespace lsmps
