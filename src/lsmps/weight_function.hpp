#pragma once

#include "core/simulation_config.hpp"

namespace lsmps {

double evaluateWeight(double distance, double support_radius, LsmpsKernelType kernel_type);

}  // namespace lsmps
