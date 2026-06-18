#pragma once

namespace lsmps {

struct CorrectionResult {
    bool applied = false;
};

class PressureCorrectionApplier {
public:
    CorrectionResult apply() const;
};

}  // namespace lsmps
