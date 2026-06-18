#include "config/config_reader.hpp"

#include "core/vector3.hpp"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lsmps {
namespace {

std::string trim(const std::string& value) {
    std::size_t begin = 0;
    while (begin < value.size() && std::isspace(static_cast<unsigned char>(value[begin]))) {
        ++begin;
    }

    std::size_t end = value.size();
    while (end > begin && std::isspace(static_cast<unsigned char>(value[end - 1]))) {
        --end;
    }

    return value.substr(begin, end - begin);
}

std::string stripComment(const std::string& line) {
    const std::size_t hash_comment = line.find('#');
    const std::size_t semicolon_comment = line.find(';');
    const std::size_t comment = std::min(hash_comment, semicolon_comment);
    if (comment == std::string::npos) {
        return line;
    }
    return line.substr(0, comment);
}

double parseDouble(const std::string& key, const std::string& value) {
    std::istringstream stream(value);
    double result = 0.0;
    stream >> result;
    if (!stream || (stream >> std::ws && !stream.eof())) {
        throw std::runtime_error("Invalid numeric value for " + key + ": " + value);
    }
    return result;
}

std::size_t parseSize(const std::string& key, const std::string& value) {
    std::istringstream stream(value);
    long long result = 0;
    stream >> result;
    if (!stream || (stream >> std::ws && !stream.eof())) {
        throw std::runtime_error("Invalid integer value for " + key + ": " + value);
    }
    if (result < 0 || static_cast<unsigned long long>(result) > std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error("Integer value out of range for " + key + ": " + value);
    }
    return static_cast<std::size_t>(result);
}

Vector3 parseVector3(const std::string& key, const std::string& value) {
    std::istringstream stream(value);
    Vector3 result;
    stream >> result.x >> result.y >> result.z;
    if (!stream || (stream >> std::ws && !stream.eof())) {
        throw std::runtime_error("Invalid Vector3 value for " + key + ": " + value);
    }
    return result;
}

std::string fullKey(const std::string& section, const std::string& key) {
    return section + "." + key;
}

bool isKnownSection(const std::string& section) {
    return section == "time" || section == "geometry" || section == "physical" ||
           section == "free_surface" || section == "linear_solver";
}

void applyTimeConfig(TimeConfig& config, const std::string& key, const std::string& value) {
    const std::string name = fullKey("time", key);
    if (key == "dt") {
        config.dt = parseDouble(name, value);
    } else if (key == "end_time") {
        config.end_time = parseDouble(name, value);
    } else if (key == "output_interval") {
        config.output_interval = parseDouble(name, value);
    } else {
        throw std::runtime_error("Unknown config key: " + name);
    }
}

void applyGeometryConfig(GeometryConfig& config, const std::string& key, const std::string& value) {
    const std::string name = fullKey("geometry", key);
    if (key == "particle_spacing") {
        config.particle_spacing = parseDouble(name, value);
    } else if (key == "support_radius") {
        config.support_radius = parseDouble(name, value);
    } else if (key == "domain_min") {
        config.domain_min = parseVector3(name, value);
    } else if (key == "domain_max") {
        config.domain_max = parseVector3(name, value);
    } else if (key == "fluid_min") {
        config.fluid_min = parseVector3(name, value);
    } else if (key == "fluid_max") {
        config.fluid_max = parseVector3(name, value);
    } else {
        throw std::runtime_error("Unknown config key: " + name);
    }
}

void applyPhysicalConfig(PhysicalConfig& config, const std::string& key, const std::string& value) {
    const std::string name = fullKey("physical", key);
    if (key == "density") {
        config.density = parseDouble(name, value);
    } else if (key == "viscosity") {
        config.viscosity = parseDouble(name, value);
    } else if (key == "gravity") {
        config.gravity = parseVector3(name, value);
    } else {
        throw std::runtime_error("Unknown config key: " + name);
    }
}

void applyFreeSurfaceConfig(FreeSurfaceConfig& config, const std::string& key, const std::string& value) {
    const std::string name = fullKey("free_surface", key);
    if (key == "neighbor_count_ratio") {
        config.neighbor_count_ratio = parseDouble(name, value);
    } else if (key == "number_density_ratio") {
        config.number_density_ratio = parseDouble(name, value);
    } else if (key == "near_surface_layers") {
        config.near_surface_layers = parseSize(name, value);
    } else {
        throw std::runtime_error("Unknown config key: " + name);
    }
}

void applyLinearSolverConfig(LinearSolverConfig& config, const std::string& key, const std::string& value) {
    const std::string name = fullKey("linear_solver", key);
    if (key == "max_iterations") {
        config.max_iterations = parseSize(name, value);
    } else if (key == "tolerance") {
        config.tolerance = parseDouble(name, value);
    } else {
        throw std::runtime_error("Unknown config key: " + name);
    }
}

void applySectionValue(
    SimulationConfig& config,
    const std::string& section,
    const std::string& key,
    const std::string& value) {
    if (section == "time") {
        applyTimeConfig(config.time, key, value);
    } else if (section == "geometry") {
        applyGeometryConfig(config.geometry, key, value);
    } else if (section == "physical") {
        applyPhysicalConfig(config.physical, key, value);
    } else if (section == "free_surface") {
        applyFreeSurfaceConfig(config.free_surface, key, value);
    } else if (section == "linear_solver") {
        applyLinearSolverConfig(config.linear_solver, key, value);
    } else {
        throw std::runtime_error("Unknown config section: " + section);
    }
}

void validateOrThrow(const SimulationConfig& config) {
    const std::vector<std::string> errors = config.validate();
    if (errors.empty()) {
        return;
    }

    std::string message = "Invalid simulation configuration:";
    for (const std::string& error : errors) {
        message += "\n  - " + error;
    }
    throw std::runtime_error(message);
}

}  // namespace

SimulationConfig ConfigReader::readFile(const std::string& path) const {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("Failed to open config file: " + path);
    }

    std::ostringstream buffer;
    buffer << file.rdbuf();
    return readString(buffer.str());
}

SimulationConfig ConfigReader::readString(const std::string& content) const {
    SimulationConfig config;
    std::string current_section;

    std::istringstream stream(content);
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(stream, line)) {
        ++line_number;
        const std::string cleaned = trim(stripComment(line));
        if (cleaned.empty()) {
            continue;
        }

        if (cleaned.front() == '[' || cleaned.back() == ']') {
            if (cleaned.front() != '[' || cleaned.back() != ']') {
                throw std::runtime_error("Invalid config section on line " + std::to_string(line_number));
            }

            current_section = trim(cleaned.substr(1, cleaned.size() - 2));
            if (current_section.empty()) {
                throw std::runtime_error("Empty config section on line " + std::to_string(line_number));
            }
            if (!isKnownSection(current_section)) {
                throw std::runtime_error(
                    "Unknown config section on line " + std::to_string(line_number) + ": " + current_section);
            }
            continue;
        }

        if (current_section.empty()) {
            throw std::runtime_error(
                "Config key outside a section on line " + std::to_string(line_number) + ": " + cleaned);
        }

        const std::size_t separator = cleaned.find('=');
        if (separator == std::string::npos) {
            throw std::runtime_error("Invalid config line " + std::to_string(line_number) + ": missing '='");
        }

        const std::string key = trim(cleaned.substr(0, separator));
        const std::string value = trim(cleaned.substr(separator + 1));
        if (key.empty() || value.empty()) {
            throw std::runtime_error("Invalid config line " + std::to_string(line_number) + ": empty key or value");
        }
        if (key.find('.') != std::string::npos) {
            throw std::runtime_error(
                "INI keys must not contain section prefixes on line " + std::to_string(line_number) + ": " + key);
        }

        try {
            applySectionValue(config, current_section, key, value);
        } catch (const std::runtime_error& error) {
            throw std::runtime_error(std::string(error.what()) + " on line " + std::to_string(line_number));
        }
    }

    validateOrThrow(config);
    return config;
}

}  // namespace lsmps
