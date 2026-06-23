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

bool parseBool(const std::string& key, const std::string& value) {
    std::string normalized = value;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });

    if (normalized == "true" || normalized == "1" || normalized == "yes" || normalized == "on") {
        return true;
    }
    if (normalized == "false" || normalized == "0" || normalized == "no" || normalized == "off") {
        return false;
    }
    throw std::runtime_error("Invalid boolean value for " + key + ": " + value);
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
    return section == "time" || section == "file" ||
           section == "geometry" || section == "physical" || section == "free_surface" ||
           section == "lsmps" || section == "linear_solver";
}

void applyTimeConfig(TimeConfig& config, const std::string& key, const std::string& value) {
    const std::string name = fullKey("time", key);
    if (key == "dt") {
        config.dt = parseDouble(name, value);
        config.initial_dt = config.dt;
    } else if (key == "start_time") {
        config.start_time = parseDouble(name, value);
    } else if (key == "end_time") {
        config.end_time = parseDouble(name, value);
    } else if (key == "initial_dt") {
        config.initial_dt = parseDouble(name, value);
        config.dt = config.initial_dt;
    } else if (key == "min_dt") {
        config.min_dt = parseDouble(name, value);
    } else if (key == "max_dt") {
        config.max_dt = parseDouble(name, value);
    } else if (key == "cfl_number") {
        config.cfl_number = parseDouble(name, value);
    } else if (key == "growth_factor") {
        config.growth_factor = parseDouble(name, value);
    } else if (key == "output_interval") {
        config.output_interval = parseDouble(name, value);
    } else {
        throw std::runtime_error("Unknown config key: " + name);
    }
}

void applyFileConfig(FileConfig& config, const std::string& key, const std::string& value) {
    const std::string name = fullKey("file", key);
    if (key == "input_directory") {
        config.input_directory = value;
    } else if (key == "input_file") {
        config.input_file = value;
    } else if (key == "fluid_particle_file") {
        config.fluid_particle_file = value;
    } else if (key == "wall_particle_file") {
        config.wall_particle_file = value;
    } else if (key == "output_directory") {
        config.output_directory = value;
    } else if (key == "output_prefix") {
        config.output_prefix = value;
    } else if (key == "write_initial_state") {
        config.write_initial_state = parseBool(name, value);
    } else if (key == "write_outputs") {
        config.write_outputs = parseBool(name, value);
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
    } else if (key == "screen_radius_factor") {
        config.screen_radius_factor = parseDouble(name, value);
    } else if (key == "wall_patch_radius_factor") {
        config.wall_patch_radius_factor = parseDouble(name, value);
    } else if (key == "particle_radius_factor") {
        config.particle_radius_factor = parseDouble(name, value);
    } else if (key == "open_threshold") {
        config.open_threshold = parseDouble(name, value);
    } else if (key == "cone_angle_degrees") {
        config.cone_angle_degrees = parseDouble(name, value);
    } else if (key == "cone_threshold") {
        config.cone_threshold = parseDouble(name, value);
    } else if (key == "min_cone_accessible_ratio") {
        config.min_cone_accessible_ratio = parseDouble(name, value);
    } else if (key == "cubed_sphere_q") {
        config.cubed_sphere_q = parseSize(name, value);
    } else if (key == "splash_max_fluid_neighbors") {
        config.splash_max_fluid_neighbors = parseSize(name, value);
    } else if (key == "splash_open_threshold") {
        config.splash_open_threshold = parseDouble(name, value);
    } else if (key == "near_surface_distance_factor") {
        config.near_surface_distance_factor = parseDouble(name, value);
    } else {
        throw std::runtime_error("Unknown config key: " + name);
    }
}

LsmpsKernelType parseKernelType(const std::string& key, const std::string& value) {
    std::string normalized = value;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(), [](unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });

    if (normalized == "linear") {
        return LsmpsKernelType::Linear;
    }
    throw std::runtime_error("Unknown kernel type for " + key + ": " + value);
}

void applyLsmpsConfig(LsmpsConfig& config, const std::string& key, const std::string& value) {
    const std::string name = fullKey("lsmps", key);
    if (key == "min_neighbors") {
        config.min_neighbors = parseSize(name, value);
    } else if (key == "eigenvalue_tolerance") {
        config.eigenvalue_tolerance = parseDouble(name, value);
    } else if (key == "condition_number_warning") {
        config.condition_number_warning = parseDouble(name, value);
    } else if (key == "condition_number_failure") {
        config.condition_number_failure = parseDouble(name, value);
    } else if (key == "diagnostics_enabled") {
        config.diagnostics_enabled = parseBool(name, value);
    } else if (key == "kernel_type") {
        config.kernel_type = parseKernelType(name, value);
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
    } else if (section == "file") {
        applyFileConfig(config.file, key, value);
    } else if (section == "geometry") {
        applyGeometryConfig(config.geometry, key, value);
    } else if (section == "physical") {
        applyPhysicalConfig(config.physical, key, value);
    } else if (section == "free_surface") {
        applyFreeSurfaceConfig(config.free_surface, key, value);
    } else if (section == "lsmps") {
        applyLsmpsConfig(config.lsmps, key, value);
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
        if (key.empty()) {
            throw std::runtime_error("Invalid config line " + std::to_string(line_number) + ": empty key");
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
