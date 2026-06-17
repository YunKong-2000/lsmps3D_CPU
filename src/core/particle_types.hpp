#pragma once

namespace lsmps {

enum class ParticleType {
    Fluid,
    Wall,
};

enum class FluidParticleState {
    Internal,
    FreeSurface,
    NearFreeSurface,
    Splash,
};

}  // namespace lsmps
