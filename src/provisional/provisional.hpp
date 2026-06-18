#pragma once

namespace lsmps {

struct ProvisionalVelocityResult {
    bool computed = false;
};

class ProvisionalVelocityCalculator {
public:
    ProvisionalVelocityResult compute() const;
};

}  // namespace lsmps
