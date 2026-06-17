#pragma once

#include <string>

namespace lsmps {

class Logger {
public:
    void info(const std::string& message) const;
    void warn(const std::string& message) const;
    void error(const std::string& message) const;
};

}  // namespace lsmps
