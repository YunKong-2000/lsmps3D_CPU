#include "io/vtk_writer.hpp"

#include "core/particle_types.hpp"

#include <cctype>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

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
    case FluidParticleState::FreeSurface:
        return 1;
    case FluidParticleState::NearFreeSurface:
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

}  // namespace

void VtkWriter::writeParticles(
    const std::string& path,
    const ParticleSet& particles,
    const std::vector<VtkScalarField>& scalar_fields,
    const std::vector<VtkVectorField>& vector_fields) const {
    const std::size_t count = particles.size();

    for (const VtkScalarField& field : scalar_fields) {
        validateFieldName(field.name);
        validateFieldSize(field.name, field.values.size(), count);
    }
    for (const VtkVectorField& field : vector_fields) {
        validateFieldName(field.name);
        validateFieldSize(field.name, field.values.size(), count);
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
    file << "POINTS " << count << " double\n";
    for (const Vector3& position : particles.positions()) {
        file << position.x << ' ' << position.y << ' ' << position.z << '\n';
    }

    file << "VERTICES " << count << ' ' << count * 2 << '\n';
    for (std::size_t index = 0; index < count; ++index) {
        file << "1 " << index << '\n';
    }

    file << "POINT_DATA " << count << '\n';

    file << "VECTORS velocity double\n";
    for (const Vector3& velocity : particles.velocities()) {
        file << velocity.x << ' ' << velocity.y << ' ' << velocity.z << '\n';
    }

    writeScalarHeader(file, "pressure", "double");
    for (const double pressure : particles.pressures()) {
        file << pressure << '\n';
    }

    writeScalarHeader(file, "particle_type", "int");
    for (const ParticleType type : particles.types()) {
        file << toVtkValue(type) << '\n';
    }

    writeScalarHeader(file, "fluid_state", "int");
    for (const FluidParticleState state : particles.fluidStates()) {
        file << toVtkValue(state) << '\n';
    }

    writeScalarHeader(file, "neighbor_count", "int");
    for (const std::size_t neighbor_count : particles.neighborCounts()) {
        file << neighbor_count << '\n';
    }

    writeScalarHeader(file, "fluid_neighbor_count", "int");
    for (const std::size_t neighbor_count : particles.fluidNeighborCounts()) {
        file << neighbor_count << '\n';
    }

    writeScalarHeader(file, "wall_neighbor_count", "int");
    for (const std::size_t neighbor_count : particles.wallNeighborCounts()) {
        file << neighbor_count << '\n';
    }

    for (const VtkScalarField& field : scalar_fields) {
        writeScalarHeader(file, field.name, "double");
        for (const double value : field.values) {
            file << value << '\n';
        }
    }

    for (const VtkVectorField& field : vector_fields) {
        file << "VECTORS " << field.name << " double\n";
        for (const Vector3& value : field.values) {
            file << value.x << ' ' << value.y << ' ' << value.z << '\n';
        }
    }
}

}  // namespace lsmps
