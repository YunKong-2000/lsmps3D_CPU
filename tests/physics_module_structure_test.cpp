#include "correction/correction.hpp"
#include "pressure_poisson/pressure_poisson.hpp"
#include "provisional/provisional.hpp"

#include <cassert>

int main() {
    const lsmps::ProvisionalVelocityCalculator provisional;
    const lsmps::PressurePoissonAssembler pressure_poisson;
    const lsmps::PressureCorrectionApplier correction;

    (void)provisional;
    (void)pressure_poisson;
    assert(!correction.apply().applied);

    return 0;
}
