#include "config/config_reader.hpp"
#include "core/simulation_config.hpp"
#include "io/logger.hpp"

#include <exception>
#include <string>

namespace {

void logConfigSummary(const lsmps::Logger& logger, const lsmps::SimulationConfig& config) {
    logger.info("dt = " + std::to_string(config.time.dt));
    logger.info("end_time = " + std::to_string(config.time.end_time));
    logger.info("output_interval = " + std::to_string(config.time.output_interval));
    logger.info("particle_spacing = " + std::to_string(config.geometry.particle_spacing));
    logger.info("support_radius = " + std::to_string(config.geometry.support_radius));
    logger.info("density = " + std::to_string(config.physical.density));
    logger.info("viscosity = " + std::to_string(config.physical.viscosity));
}

}  // namespace

int main(int argc, char** argv) {
    lsmps::Logger logger;
    logger.info("LSMPS 3D CPU solver skeleton");
    logger.info("Build system and module layout are ready.");

    try {
        lsmps::SimulationConfig config;
        if (argc > 2) {
            logger.error("Usage: lsmps3d [config_file]");
            return 1;
        }
        if (argc == 2) {
            const lsmps::ConfigReader reader;
            config = reader.readFile(argv[1]);
            logger.info(std::string("Loaded config file: ") + argv[1]);
        } else {
            logger.info("Using default simulation configuration.");
        }
        logConfigSummary(logger, config);
    } catch (const std::exception& error) {
        logger.error(error.what());
        return 1;
    }

    return 0;
}
