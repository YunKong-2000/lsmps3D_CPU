#include "lsmps/weight_function.hpp"

#include <algorithm>

namespace lsmps {

double evaluateWeight(double distance, double support_radius, LsmpsKernelType kernel_type) {
    if (distance <= 0.0 || distance > support_radius) {
        return 0.0;
    }

    switch (kernel_type) {
    case LsmpsKernelType::Linear:
        return std::max(0.0, pow(1.0 - distance / support_radius, 2));
    }

    return 0.0;
}

}  // namespace lsmps
