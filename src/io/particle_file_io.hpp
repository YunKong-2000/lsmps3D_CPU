#pragma once

#include "core/particle_set.hpp"

#include <string>

namespace lsmps {

class ParticleFileReader {
public:
    void readFluidParticles(const std::string& path, ParticleSet& particles) const;
    void readWallParticles(const std::string& path, ParticleSet& particles) const;
};

class ParticleFileWriter {
public:
    void writeFluidParticles(const std::string& path, const ParticleSet& particles) const;
    void writeWallParticles(const std::string& path, const ParticleSet& particles) const;
};

}  // namespace lsmps
