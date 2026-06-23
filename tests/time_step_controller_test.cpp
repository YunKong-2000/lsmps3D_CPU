#include "time_integration/time_step_controller.hpp"

#include <cassert>
#include <cmath>

int main() {
    lsmps::TimeStepControlConfig config;
    config.start_time = 0.0;
    config.end_time = 1.0;
    config.initial_dt = 0.01;
    config.min_dt = 0.001;
    config.max_dt = 0.08;
    config.cfl_number = 0.2;
    config.growth_factor = 2.0;
    config.output_interval = 0.05;

    lsmps::TimeStepController controller(config, 0.1);

    const lsmps::TimeStepDecision grown = controller.chooseDt(0.0, 0.1);
    assert(std::abs(grown.dt - 0.02) < 1.0e-14);
    assert(std::abs(grown.cfl_limited_dt - 0.2) < 1.0e-14);
    controller.updateAfterStep(grown);

    const lsmps::TimeStepDecision output_limited = controller.chooseDt(0.02, 0.1);
    assert(std::abs(output_limited.dt - 0.03) < 1.0e-14);
    assert(output_limited.limited_by_output_time);
    controller.updateAfterStep(output_limited);
    assert(controller.shouldWriteOutput(0.05));
    controller.advanceOutputTime(0.05);
    assert(!controller.shouldWriteOutput(0.05));

    const lsmps::TimeStepDecision cfl_limited = controller.chooseDt(0.05, 10.0);
    assert(std::abs(cfl_limited.cfl_limited_dt - 0.002) < 1.0e-14);
    assert(std::abs(cfl_limited.dt - 0.002) < 1.0e-14);
    assert(cfl_limited.limited_by_cfl);

    const lsmps::TimeStepDecision min_limited = controller.chooseDt(0.052, 1000.0);
    assert(std::abs(min_limited.dt - config.min_dt) < 1.0e-14);
    assert(min_limited.limited_by_min_dt);

    const lsmps::TimeStepDecision end_limited = controller.chooseDt(0.9995, 0.0);
    assert(std::abs(end_limited.dt - 0.0005) < 1.0e-14);
    assert(end_limited.limited_by_end_time);

    return 0;
}
