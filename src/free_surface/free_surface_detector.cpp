#include "free_surface/free_surface_detector.hpp"

#include "core/particle_types.hpp"
#include "core/vector3.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace lsmps {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

Vector3 normalize(const Vector3& value) {
    const double length = norm(value);
    if (length <= std::numeric_limits<double>::epsilon()) {
        return {};
    }
    return value / length;
}

Vector3 cubedSphereDirection(int face, double u, double v) {
    switch (face) {
    case 0:
        return normalize({1.0, u, v});
    case 1:
        return normalize({-1.0, u, v});
    case 2:
        return normalize({u, 1.0, v});
    case 3:
        return normalize({u, -1.0, v});
    case 4:
        return normalize({u, v, 1.0});
    default:
        return normalize({u, v, -1.0});
    }
}

bool validWallNormal(const Vector3& normal) {
    return normSquared(normal) > 0.0;
}

}  // namespace

FreeSurfaceDetector::FreeSurfaceDetector(const FreeSurfaceConfig& config, double particle_spacing)
    : config_(config),
      particle_spacing_(particle_spacing),
      screen_radius_(config.screen_radius_factor * particle_spacing),
      wall_patch_radius_(config.wall_patch_radius_factor * particle_spacing),
      particle_radius_(config.particle_radius_factor * particle_spacing),
      cone_cosine_(std::cos(config.cone_angle_degrees * pi / 180.0)) {
    if (particle_spacing_ <= 0.0) {
        throw std::runtime_error("FreeSurfaceDetector particle spacing must be positive");
    }
    if (config_.cubed_sphere_q == 0) {
        throw std::runtime_error("FreeSurfaceDetector cubed_sphere_q must be positive");
    }
    buildCubedSphere(config_.cubed_sphere_q);
    buildConeTable();
}

FreeSurfaceDiagnostics FreeSurfaceDetector::detect(
    ParticleSet& particles,
    const TypedNeighborList& neighbors) const {
    if (neighbors.fluid.size() != particles.size() || neighbors.wall.size() != particles.size()) {
        throw std::runtime_error("Free-surface neighbor list size must match particle count");
    }

    const std::size_t particle_count = particles.size();
    const std::size_t direction_count = directions_.size();
    FreeSurfaceDiagnostics diagnostics;
    diagnostics.open_ratio.assign(particle_count, 0.0);
    diagnostics.cone_ratio.assign(particle_count, 0.0);
    diagnostics.accessible_area_ratio.assign(particle_count, 0.0);
    diagnostics.reason_code.assign(particle_count, 0);

    std::vector<FluidParticleState> primary_states(particle_count, FluidParticleState::Internal);
    const double total_area = std::max(
        std::numeric_limits<double>::epsilon(),
        [&]() {
            double area = 0.0;
            for (const double weight : direction_weights_) {
                area += weight;
            }
            return area;
        }());

    for (std::size_t i = 0; i < particle_count; ++i) {
        if (particles.types()[i] != ParticleType::Fluid) {
            continue;
        }

        std::vector<unsigned char> solid(direction_count, 0);
        std::vector<unsigned char> shadow(direction_count, 0);
        const Vector3& origin = particles.positions()[i];

        for (const std::size_t wall_index : neighbors.wall[i]) {
            const Vector3 wall_normal = normalize(particles.wallNormals()[wall_index]);
            if (!validWallNormal(wall_normal)) {
                continue;
            }
            const Vector3 wall_offset = particles.positions()[wall_index] - origin;

            for (std::size_t m = 0; m < direction_count; ++m) {
                if (solid[m]) {
                    continue;
                }
                const Vector3& direction = directions_[m];
                const double den = dot(direction, wall_normal);
                if (den >= 0.0) {
                    continue;
                }

                const double t = dot(wall_offset, wall_normal) / den;
                if (t <= 0.0 || t >= screen_radius_) {
                    continue;
                }

                const Vector3 hit = origin + direction * t;
                const Vector3 hit_offset = hit - particles.positions()[wall_index];
                const Vector3 tangent = hit_offset - wall_normal * dot(hit_offset, wall_normal);
                if (norm(tangent) <= wall_patch_radius_) {
                    solid[m] = 1;
                }
            }
        }

        for (const std::size_t fluid_index : neighbors.fluid[i]) {
            const Vector3 offset = particles.positions()[fluid_index] - origin;
            const double distance = norm(offset);
            if (distance <= std::numeric_limits<double>::epsilon()) {
                continue;
            }
            const Vector3 direction_to_neighbor = offset / distance;
            const double radius_ratio = std::min(1.0, particle_radius_ / distance);
            const double cos_beta = std::sqrt(std::max(0.0, 1.0 - radius_ratio * radius_ratio));

            for (std::size_t m = 0; m < direction_count; ++m) {
                if (!solid[m] && dot(directions_[m], direction_to_neighbor) >= cos_beta) {
                    shadow[m] = 1;
                }
            }
        }

        double accessible_area = 0.0;
        double open_area = 0.0;
        for (std::size_t m = 0; m < direction_count; ++m) {
            if (!solid[m]) {
                accessible_area += direction_weights_[m];
                if (!shadow[m]) {
                    open_area += direction_weights_[m];
                }
            }
        }

        if (accessible_area <= std::numeric_limits<double>::epsilon()) {
            diagnostics.reason_code[i] = 5;
            continue;
        }

        const double open_ratio = open_area / accessible_area;
        double cone_ratio = 0.0;
        for (std::size_t cone = 0; cone < cone_direction_indices_.size(); ++cone) {
            double cone_accessible_area = 0.0;
            double cone_open_area = 0.0;
            for (const std::size_t m : cone_direction_indices_[cone]) {
                if (!solid[m]) {
                    cone_accessible_area += direction_weights_[m];
                    if (!shadow[m]) {
                        cone_open_area += direction_weights_[m];
                    }
                }
            }

            if (cone_areas_[cone] <= std::numeric_limits<double>::epsilon()) {
                continue;
            }
            const double accessible_ratio = cone_accessible_area / cone_areas_[cone];
            if (accessible_ratio < config_.min_cone_accessible_ratio ||
                cone_accessible_area <= std::numeric_limits<double>::epsilon()) {
                continue;
            }
            cone_ratio = std::max(cone_ratio, cone_open_area / cone_accessible_area);
        }

        diagnostics.open_ratio[i] = open_ratio;
        diagnostics.cone_ratio[i] = cone_ratio;
        diagnostics.accessible_area_ratio[i] = accessible_area / total_area;

        if (neighbors.fluid[i].size() <= config_.splash_max_fluid_neighbors && neighbors.wall[i].empty() &&
            open_ratio >= config_.splash_open_threshold) {
            primary_states[i] = FluidParticleState::Splash;
            diagnostics.reason_code[i] = 4;
        } else if (open_ratio >= config_.open_threshold && cone_ratio >= config_.cone_threshold) {
            primary_states[i] = FluidParticleState::FreeSurface;
            diagnostics.reason_code[i] = 1;
        }
    }

    std::vector<FluidParticleState>& states = particles.fluidStates();
    const double near_surface_distance = config_.near_surface_distance_factor * particle_spacing_;
    const double near_surface_distance_squared = near_surface_distance * near_surface_distance;
    for (std::size_t i = 0; i < particle_count; ++i) {
        if (particles.types()[i] != ParticleType::Fluid) {
            continue;
        }
        if (primary_states[i] == FluidParticleState::FreeSurface || primary_states[i] == FluidParticleState::Splash) {
            states[i] = primary_states[i];
            continue;
        }

        bool near_surface = false;
        for (const std::size_t neighbor : neighbors.fluid[i]) {
            if (primary_states[neighbor] == FluidParticleState::FreeSurface) {
                const Vector3 offset = particles.positions()[neighbor] - particles.positions()[i];
                if (normSquared(offset) < near_surface_distance_squared) {
                    near_surface = true;
                    break;
                }
            }
        }
        states[i] = near_surface ? FluidParticleState::NearFreeSurface : FluidParticleState::Internal;
        if (near_surface) {
            diagnostics.reason_code[i] = 3;
        }
    }

    return diagnostics;
}

const std::vector<Vector3>& FreeSurfaceDetector::directions() const noexcept {
    return directions_;
}

const std::vector<double>& FreeSurfaceDetector::directionWeights() const noexcept {
    return direction_weights_;
}

void FreeSurfaceDetector::buildCubedSphere(std::size_t q) {
    directions_.clear();
    direction_weights_.clear();
    directions_.reserve(6 * q * q);
    direction_weights_.reserve(6 * q * q);

    const double du = 2.0 / static_cast<double>(q);
    for (int face = 0; face < 6; ++face) {
        for (std::size_t iy = 0; iy < q; ++iy) {
            for (std::size_t ix = 0; ix < q; ++ix) {
                const double u = -1.0 + (static_cast<double>(ix) + 0.5) * du;
                const double v = -1.0 + (static_cast<double>(iy) + 0.5) * du;
                directions_.push_back(cubedSphereDirection(face, u, v));
                const double length = std::sqrt(1.0 + u * u + v * v);
                direction_weights_.push_back((du * du) / (length * length * length));
            }
        }
    }
}

void FreeSurfaceDetector::buildConeTable() {
    cone_direction_indices_.assign(directions_.size(), {});
    cone_areas_.assign(directions_.size(), 0.0);

    for (std::size_t q = 0; q < directions_.size(); ++q) {
        for (std::size_t m = 0; m < directions_.size(); ++m) {
            if (dot(directions_[q], directions_[m]) >= cone_cosine_) {
                cone_direction_indices_[q].push_back(m);
                cone_areas_[q] += direction_weights_[m];
            }
        }
    }
}

}  // namespace lsmps
