#include "logger.hpp"

#include <iostream>

namespace lsmps {

void Logger::info(const std::string& message) const {
    std::cout << "[info] " << message << '\n';
}

void Logger::warn(const std::string& message) const {
    std::cout << "[warn] " << message << '\n';
}

void Logger::error(const std::string& message) const {
    std::cerr << "[error] " << message << '\n';
}

}  // namespace lsmps
