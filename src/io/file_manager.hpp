#pragma once

#include "core/simulation_config.hpp"

#include <cstddef>
#include <string>

namespace lsmps {

class FileManager {
public:
    explicit FileManager(FileConfig config = {});

    const FileConfig& config() const noexcept;

    std::string inputPath() const;
    std::string fluidParticlePath() const;
    std::string wallParticlePath() const;
    std::string outputPath(const std::string& tag) const;
    std::string initialOutputPath() const;
    std::string stepOutputPath(std::size_t output_index, std::size_t step) const;
    std::string fluidOutputPath(std::size_t output_index) const;
    std::string wallOutputPath(std::size_t output_index) const;

private:
    FileConfig config_;
};

}  // namespace lsmps
