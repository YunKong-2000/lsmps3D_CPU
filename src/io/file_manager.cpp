#include "io/file_manager.hpp"

#include <filesystem>
#include <iomanip>
#include <sstream>
#include <utility>

namespace lsmps {
namespace {

std::string formatStepTag(std::size_t output_index, std::size_t step) {
    std::ostringstream stream;
    stream << std::setfill('0') << std::setw(5) << output_index << "_step_" << step;
    return stream.str();
}

}  // namespace

FileManager::FileManager(FileConfig config)
    : config_(std::move(config)) {}

const FileConfig& FileManager::config() const noexcept {
    return config_;
}

std::string FileManager::inputPath() const {
    if (config_.input_file.empty()) {
        return {};
    }
    const std::filesystem::path file(config_.input_file);
    if (file.is_absolute() || config_.input_directory.empty()) {
        return file.string();
    }
    return (std::filesystem::path(config_.input_directory) / file).string();
}

std::string FileManager::outputPath(const std::string& tag) const {
    const std::filesystem::path directory(config_.output_directory);
    const std::filesystem::path file = config_.output_prefix + "_" + tag + ".vtk";
    return (directory / file).string();
}

std::string FileManager::initialOutputPath() const {
    return outputPath("initial");
}

std::string FileManager::stepOutputPath(std::size_t output_index, std::size_t step) const {
    return outputPath(formatStepTag(output_index, step));
}

}  // namespace lsmps
