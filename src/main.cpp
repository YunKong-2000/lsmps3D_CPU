#include "config/config_reader.hpp"
#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "core/vector3.hpp"
#include "io/logger.hpp"
#include "time_integration/time_stepper.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace {

int gridPointCount(double min_value, double max_value, double spacing) {
    return static_cast<int>(std::round((max_value - min_value) / spacing)) + 1;
}

double coordinate(double min_value, int index, double spacing) {
    return min_value + static_cast<double>(index) * spacing;
}

void logConfigSummary(const lsmps::Logger& logger, const lsmps::SimulationConfig& config) {
    logger.info("dt = " + std::to_string(config.time.dt));
    logger.info("end_time = " + std::to_string(config.time.end_time));
    logger.info("output_interval = " + std::to_string(config.time.output_interval));
    logger.info("particle_spacing = " + std::to_string(config.geometry.particle_spacing));
    logger.info("support_radius = " + std::to_string(config.geometry.support_radius));
    logger.info("density = " + std::to_string(config.physical.density));
    logger.info("viscosity = " + std::to_string(config.physical.viscosity));
}

lsmps::SimulationConfig compactDemoConfig() {
    lsmps::SimulationConfig config;
    config.time.dt = 0.002;
    config.time.end_time = 0.02;
    config.time.output_interval = 0.01;
    config.geometry.particle_spacing = 0.2;
    config.geometry.support_radius = 3.1 * config.geometry.particle_spacing;
    config.geometry.domain_min = {0.0, 0.0, 0.0};
    config.geometry.domain_max = {1.0, 1.0, 1.0};
    config.geometry.fluid_min = {0.0, 0.0, 0.0};
    config.geometry.fluid_max = {1.0, 0.6, 1.0};
    config.physical.density = 1000.0;
    config.physical.viscosity = 0.0;
    config.physical.gravity = {0.0, -9.81, 0.0};
    config.lsmps.condition_number_warning = 1.0e10;
    config.lsmps.condition_number_failure = 1.0e14;
    config.linear_solver.max_iterations = 1000;
    config.linear_solver.tolerance = 1.0e-10;
    return config;
}

void addFluidParticles(lsmps::ParticleSet& particles, const lsmps::SimulationConfig& config) {
    const double spacing = config.geometry.particle_spacing;
    const int nx = gridPointCount(config.geometry.fluid_min.x, config.geometry.fluid_max.x, spacing);
    const int ny = gridPointCount(config.geometry.fluid_min.y, config.geometry.fluid_max.y, spacing);
    const int nz = gridPointCount(config.geometry.fluid_min.z, config.geometry.fluid_max.z, spacing);

    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                particles.addFluidParticle(
                    {
                        coordinate(config.geometry.fluid_min.x, x, spacing),
                        coordinate(config.geometry.fluid_min.y, y, spacing),
                        coordinate(config.geometry.fluid_min.z, z, spacing),
                    },
                    {},
                    lsmps::FluidParticleState::Internal,
                    0.0,
                    config.physical.density);
            }
        }
    }
}

void addTankWallParticles(lsmps::ParticleSet& particles, const lsmps::SimulationConfig& config) {
    const double spacing = config.geometry.particle_spacing;
    const int nx = gridPointCount(config.geometry.domain_min.x, config.geometry.domain_max.x, spacing);
    const int ny = gridPointCount(config.geometry.domain_min.y, config.geometry.domain_max.y, spacing);
    const int nz = gridPointCount(config.geometry.domain_min.z, config.geometry.domain_max.z, spacing);

    for (int z = -1; z <= nz; ++z) {
        for (int y = -1; y < ny; ++y) {
            for (int x = -1; x <= nx; ++x) {
                lsmps::Vector3 normal;
                bool is_wall = false;
                if (x == -1) {
                    normal += {1.0, 0.0, 0.0};
                    is_wall = true;
                }
                if (x == nx) {
                    normal += {-1.0, 0.0, 0.0};
                    is_wall = true;
                }
                if (y == -1) {
                    normal += {0.0, 1.0, 0.0};
                    is_wall = true;
                }
                if (z == -1) {
                    normal += {0.0, 0.0, 1.0};
                    is_wall = true;
                }
                if (z == nz) {
                    normal += {0.0, 0.0, -1.0};
                    is_wall = true;
                }

                if (!is_wall) {
                    continue;
                }

                const double length = lsmps::norm(normal);
                if (length > 0.0) {
                    normal /= length;
                }
                particles.addWallParticle(
                    {
                        coordinate(config.geometry.domain_min.x, x, spacing),
                        coordinate(config.geometry.domain_min.y, y, spacing),
                        coordinate(config.geometry.domain_min.z, z, spacing),
                    },
                    {},
                    0.0,
                    config.physical.density,
                    normal);
            }
        }
    }
}

lsmps::ParticleSet createHydrostaticBoxParticles(const lsmps::SimulationConfig& config) {
    lsmps::ParticleSet particles;
    const double spacing = config.geometry.particle_spacing;
    const int fluid_nx = gridPointCount(config.geometry.fluid_min.x, config.geometry.fluid_max.x, spacing);
    const int fluid_ny = gridPointCount(config.geometry.fluid_min.y, config.geometry.fluid_max.y, spacing);
    const int fluid_nz = gridPointCount(config.geometry.fluid_min.z, config.geometry.fluid_max.z, spacing);
    const int domain_nx = gridPointCount(config.geometry.domain_min.x, config.geometry.domain_max.x, spacing);
    const int domain_ny = gridPointCount(config.geometry.domain_min.y, config.geometry.domain_max.y, spacing);
    const int domain_nz = gridPointCount(config.geometry.domain_min.z, config.geometry.domain_max.z, spacing);
    particles.reserve(
        static_cast<std::size_t>(fluid_nx * fluid_ny * fluid_nz) +
        static_cast<std::size_t>((domain_nx + 2) * (domain_ny + 1) * (domain_nz + 2)));
    addFluidParticles(particles, config);
    addTankWallParticles(particles, config);
    return particles;
}

void logDiagnostics(const lsmps::Logger& logger, const std::vector<lsmps::TimeStepDiagnostics>& history) {
    for (const lsmps::TimeStepDiagnostics& diagnostics : history) {
        std::ostringstream message;
        message << "step=" << diagnostics.step
                << " time=" << std::setprecision(6) << diagnostics.time
                << " max_velocity=" << diagnostics.max_velocity
                << " cfl=" << diagnostics.cfl_number
                << " ppe_iterations=" << diagnostics.pressure_iterations
                << " ppe_residual=" << diagnostics.pressure_residual
                << " free_surface=" << diagnostics.free_surface_count
                << " near_surface=" << diagnostics.near_free_surface_count;
        logger.info(message.str());
    }
}

}  // namespace

int main(int argc, char** argv) {
    lsmps::Logger logger;
    logger.info("LSMPS 3D CPU solver");

    try {
        lsmps::SimulationConfig config = compactDemoConfig();
        if (argc > 2) {
            logger.error("Usage: lsmps3d [config_file]");
            return 1;
        }
        if (argc == 2) {
            const lsmps::ConfigReader reader;
            config = reader.readFile(argv[1]);
            logger.info(std::string("Loaded config file: ") + argv[1]);
        } else {
            logger.info("Using compact built-in hydrostatic demonstration configuration.");
        }
        logConfigSummary(logger, config);

        lsmps::ParticleSet particles = createHydrostaticBoxParticles(config);
        logger.info("Created hydrostatic box particles: " + std::to_string(particles.size()));

        lsmps::TimeStepperOptions options;
        options.output_directory = "output/main_hydrostatic";
        options.output_prefix = "hydrostatic";

        lsmps::TimeStepper stepper(config, options);
        const std::vector<lsmps::TimeStepDiagnostics> history = stepper.run(particles);
        logDiagnostics(logger, history);
        logger.info("Completed time integration. VTK output directory: " + options.output_directory);
    } catch (const std::exception& error) {
        logger.error(error.what());
        return 1;
    }

    return 0;
}
