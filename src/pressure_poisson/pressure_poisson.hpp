#pragma once

namespace lsmps {

struct PressurePoissonResult {
    bool solved = false;
};

class PressurePoissonAssembler {
public:
    PressurePoissonResult assembleAndSolve() const;
};

}  // namespace lsmps
