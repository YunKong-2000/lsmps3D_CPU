#include "io/particle_file_io.hpp"

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace lsmps {
namespace {

void ensureParentDirectoryExists(const std::string& path) {
    const std::filesystem::path output_path(path);
    const std::filesystem::path parent = output_path.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

void skipHeader(std::istream& input, const std::string& path, std::size_t& count) {
    std::string header;
    if (!std::getline(input, header)) {
        throw std::runtime_error("Failed to read particle file header: " + path);
    }
    if (!(input >> count)) {
        throw std::runtime_error("Failed to read particle count: " + path);
    }
}

}  // namespace

void ParticleFileReader::readFluidParticles(const std::string& path, ParticleSet& particles) const {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open fluid particle file: " + path);
    }

    std::size_t count = 0;
    skipHeader(input, path, count);
    for (std::size_t index = 0; index < count; ++index) {
        Vector3 position;
        Vector3 velocity;
        double pressure = 0.0;
        double density = 0.0;
        if (!(input >> position.x >> position.y >> position.z
              >> velocity.x >> velocity.y >> velocity.z
              >> pressure >> density)) {
            throw std::runtime_error("Failed to read fluid particle entry from: " + path);
        }
        particles.addFluidParticle(position, velocity, FluidParticleState::Internal, pressure, density);
    }
}

void ParticleFileReader::readWallParticles(const std::string& path, ParticleSet& particles) const {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open wall particle file: " + path);
    }

    std::size_t count = 0;
    skipHeader(input, path, count);
    for (std::size_t index = 0; index < count; ++index) {
        Vector3 position;
        Vector3 velocity;
        Vector3 normal;
        double pressure = 0.0;
        double density = 0.0;
        if (!(input >> position.x >> position.y >> position.z
              >> velocity.x >> velocity.y >> velocity.z
              >> pressure >> density
              >> normal.x >> normal.y >> normal.z)) {
            throw std::runtime_error("Failed to read wall particle entry from: " + path);
        }
        particles.addWallParticle(position, velocity, pressure, density, normal);
    }
}

void ParticleFileWriter::writeFluidParticles(const std::string& path, const ParticleSet& particles) const {
    ensureParentDirectoryExists(path);
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("Failed to open fluid particle output file: " + path);
    }

    std::size_t count = 0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isFluid(i)) {
            ++count;
        }
    }

    output << std::setprecision(17);
    output << "# x y z vx vy vz pressure density\n";
    output << count << '\n';
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!particles.isFluid(i)) {
            continue;
        }
        const Vector3& x = particles.positions()[i];
        const Vector3& v = particles.velocities()[i];
        output << x.x << ' ' << x.y << ' ' << x.z << ' '
               << v.x << ' ' << v.y << ' ' << v.z << ' '
               << particles.pressures()[i] << ' ' << particles.densities()[i] << '\n';
    }
}

void ParticleFileWriter::writeWallParticles(const std::string& path, const ParticleSet& particles) const {
    ensureParentDirectoryExists(path);
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error("Failed to open wall particle output file: " + path);
    }

    std::size_t count = 0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (particles.isWall(i)) {
            ++count;
        }
    }

    output << std::setprecision(17);
    output << "# x y z vx vy vz pressure density nx ny nz\n";
    output << count << '\n';
    for (std::size_t i = 0; i < particles.size(); ++i) {
        if (!particles.isWall(i)) {
            continue;
        }
        const Vector3& x = particles.positions()[i];
        const Vector3& v = particles.velocities()[i];
        const Vector3& n = particles.wallNormals()[i];
        output << x.x << ' ' << x.y << ' ' << x.z << ' '
               << v.x << ' ' << v.y << ' ' << v.z << ' '
               << particles.pressures()[i] << ' ' << particles.densities()[i] << ' '
               << n.x << ' ' << n.y << ' ' << n.z << '\n';
    }
}

}  // namespace lsmps
