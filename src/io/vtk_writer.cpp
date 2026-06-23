#include "io/vtk_writer.hpp"

#include "core/particle_types.hpp"

#include <cctype>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace lsmps {
namespace {

int toVtkValue(ParticleType type) {
    switch (type) {
    case ParticleType::Fluid:
        return 0;
    case ParticleType::Wall:
        return 1;
    }
    return -1;
}

int toVtkValue(FluidParticleState state) {
    switch (state) {
    case FluidParticleState::Internal:
        return 0;
    case FluidParticleState::NearFreeSurface:
        return 1;
    case FluidParticleState::FreeSurface:
        return 2;
    case FluidParticleState::Splash:
        return 3;
    }
    return -1;
}

void validateFieldName(const std::string& name) {
    if (name.empty()) {
        throw std::runtime_error("VTK field name must not be empty");
    }
    for (const char value : name) {
        if (std::isspace(static_cast<unsigned char>(value))) {
            throw std::runtime_error("VTK field name must not contain whitespace: " + name);
        }
    }
}

void validateFieldSize(const std::string& name, std::size_t field_size, std::size_t particle_count) {
    if (field_size != particle_count) {
        throw std::runtime_error(
            "VTK field size mismatch for " + name + ": expected " + std::to_string(particle_count) +
            ", got " + std::to_string(field_size));
    }
}

void ensureParentDirectoryExists(const std::string& path) {
    const std::filesystem::path output_path(path);
    const std::filesystem::path parent = output_path.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

void writeScalarHeader(std::ofstream& file, const std::string& name, const std::string& type) {
    file << "SCALARS " << name << ' ' << type << " 1\n";
    file << "LOOKUP_TABLE default\n";
}

std::vector<std::size_t> allParticleIndices(const ParticleSet& particles) {
    std::vector<std::size_t> indices(particles.size());
    for (std::size_t i = 0; i < particles.size(); ++i) {
        indices[i] = i;
    }
    return indices;
}

std::vector<std::size_t> particleIndicesByType(const ParticleSet& particles, ParticleType type) {
    std::vector<std::size_t> indices;
    indices.reserve(particles.size());
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.types()[i] == type) {
            indices.push_back(i);
        }
    }
    return indices;
}

void writeParticleSubset(
    const std::string& path,
    const ParticleSet& particles,
    const std::vector<std::size_t>& indices,
    const VtkBuiltInFieldOptions& built_in_fields,
    const std::vector<VtkScalarField>& scalar_fields,
    const std::vector<VtkVectorField>& vector_fields) {
    const std::size_t particle_count = particles.size();
    for (const VtkScalarField& field : scalar_fields) {
        validateFieldName(field.name);
        validateFieldSize(field.name, field.values.size(), particle_count);
    }
    for (const VtkVectorField& field : vector_fields) {
        validateFieldName(field.name);
        validateFieldSize(field.name, field.values.size(), particle_count);
    }

    ensureParentDirectoryExists(path);
    std::ofstream file(path);
    if (!file) {
        throw std::runtime_error("Failed to open VTK output file: " + path);
    }

    file << std::setprecision(17);
    file << "# vtk DataFile Version 3.0\n";
    file << "LSMPS 3D particle data\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";
    file << "POINTS " << indices.size() << " double\n";
    for (const std::size_t index : indices) {
        const Vector3& position = particles.positions()[index];
        file << position.x << ' ' << position.y << ' ' << position.z << '\n';
    }

    file << "VERTICES " << indices.size() << ' ' << indices.size() * 2 << '\n';
    for (std::size_t index = 0; index < indices.size(); ++index) {
        file << "1 " << index << '\n';
    }

    file << "POINT_DATA " << indices.size() << '\n';

    if (built_in_fields.velocity) {
        file << "VECTORS velocity double\n";
        for (const std::size_t index : indices) {
            const Vector3& velocity = particles.velocities()[index];
            file << velocity.x << ' ' << velocity.y << ' ' << velocity.z << '\n';
        }
    }

    if (built_in_fields.pressure) {
        writeScalarHeader(file, "pressure", "double");
        for (const std::size_t index : indices) {
            file << particles.pressures()[index] << '\n';
        }
    }

    if (built_in_fields.particle_type) {
        writeScalarHeader(file, "particle_type", "int");
        for (const std::size_t index : indices) {
            file << toVtkValue(particles.types()[index]) << '\n';
        }
    }

    if (built_in_fields.fluid_state) {
        writeScalarHeader(file, "fluid_state", "int");
        for (const std::size_t index : indices) {
            file << toVtkValue(particles.fluidStates()[index]) << '\n';
        }
    }

    if (built_in_fields.neighbor_count) {
        writeScalarHeader(file, "neighbor_count", "int");
        for (const std::size_t index : indices) {
            file << particles.neighborCounts()[index] << '\n';
        }
    }

    if (built_in_fields.fluid_neighbor_count) {
        writeScalarHeader(file, "fluid_neighbor_count", "int");
        for (const std::size_t index : indices) {
            file << particles.fluidNeighborCounts()[index] << '\n';
        }
    }

    if (built_in_fields.wall_neighbor_count) {
        writeScalarHeader(file, "wall_neighbor_count", "int");
        for (const std::size_t index : indices) {
            file << particles.wallNeighborCounts()[index] << '\n';
        }
    }

    for (const VtkScalarField& field : scalar_fields) {
        writeScalarHeader(file, field.name, "double");
        for (const std::size_t index : indices) {
            file << field.values[index] << '\n';
        }
    }

    for (const VtkVectorField& field : vector_fields) {
        file << "VECTORS " << field.name << " double\n";
        for (const std::size_t index : indices) {
            const Vector3& value = field.values[index];
            file << value.x << ' ' << value.y << ' ' << value.z << '\n';
        }
    }
}

}  // namespace

void VtkWriter::writeParticles(
    const std::string& path,
    const ParticleSet& particles,
    const std::vector<VtkScalarField>& scalar_fields,
    const std::vector<VtkVectorField>& vector_fields) const {
    writeParticleSubset(path, particles, allParticleIndices(particles), {}, scalar_fields, vector_fields);
}

void VtkWriter::writeParticlesByType(
    const std::string& path,
    const ParticleSet& particles,
    ParticleType type,
    const std::vector<VtkScalarField>& scalar_fields,
    const std::vector<VtkVectorField>& vector_fields) const {
    writeParticleSubset(path, particles, particleIndicesByType(particles, type), {}, scalar_fields, vector_fields);
}

void VtkWriter::writeParticlesByType(
    const std::string& path,
    const ParticleSet& particles,
    ParticleType type,
    const VtkBuiltInFieldOptions& built_in_fields,
    const std::vector<VtkScalarField>& scalar_fields,
    const std::vector<VtkVectorField>& vector_fields) const {
    writeParticleSubset(
        path,
        particles,
        particleIndicesByType(particles, type),
        built_in_fields,
        scalar_fields,
        vector_fields);
}

}  // namespace lsmps
