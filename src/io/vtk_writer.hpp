#pragma once

#include "core/particle_set.hpp"
#include "core/vector3.hpp"

#include <string>
#include <vector>

namespace lsmps {

struct VtkScalarField {
    std::string name;
    std::vector<double> values;
};

struct VtkVectorField {
    std::string name;
    std::vector<Vector3> values;
};

class VtkWriter {
public:
    void writeParticles(
        const std::string& path,
        const ParticleSet& particles,
        const std::vector<VtkScalarField>& scalar_fields = {},
        const std::vector<VtkVectorField>& vector_fields = {}) const;
};

}  // namespace lsmps
