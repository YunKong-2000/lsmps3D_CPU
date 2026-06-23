#pragma once

#include "core/particle_set.hpp"
#include "core/particle_types.hpp"
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

struct VtkBuiltInFieldOptions {
    bool velocity = true;
    bool pressure = true;
    bool particle_type = true;
    bool fluid_state = true;
    bool neighbor_count = true;
    bool fluid_neighbor_count = true;
    bool wall_neighbor_count = true;
};

class VtkWriter {
public:
    void writeParticles(
        const std::string& path,
        const ParticleSet& particles,
        const std::vector<VtkScalarField>& scalar_fields = {},
        const std::vector<VtkVectorField>& vector_fields = {}) const;

    void writeParticlesByType(
        const std::string& path,
        const ParticleSet& particles,
        ParticleType type,
        const std::vector<VtkScalarField>& scalar_fields = {},
        const std::vector<VtkVectorField>& vector_fields = {}) const;

    void writeParticlesByType(
        const std::string& path,
        const ParticleSet& particles,
        ParticleType type,
        const VtkBuiltInFieldOptions& built_in_fields,
        const std::vector<VtkScalarField>& scalar_fields = {},
        const std::vector<VtkVectorField>& vector_fields = {}) const;
};

}  // namespace lsmps
