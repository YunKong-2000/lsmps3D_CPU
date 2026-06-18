#pragma once

namespace lsmps {

enum class ParticleType {
    Fluid,
    Wall,
};

enum class FluidParticleState {
    Internal = 0,
    NearFreeSurface = 1,
    FreeSurface = 2,
    Splash = 3,
};

}  // namespace lsmps
