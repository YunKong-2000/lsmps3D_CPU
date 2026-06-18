#pragma once

#include "core/simulation_config.hpp"

#include <string>

namespace lsmps {

class ConfigReader {
public:
    SimulationConfig readFile(const std::string& path) const;
    SimulationConfig readString(const std::string& content) const;
};

}  // namespace lsmps
