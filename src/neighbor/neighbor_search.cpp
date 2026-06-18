#include "neighbor/neighbor_search.hpp"

#include "core/particle_types.hpp"
#include "core/vector3.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>

namespace lsmps {
namespace {

void validateSupportRadius(double support_radius) {
    if (support_radius <= 0.0) {
        throw std::runtime_error("NeighborSearch support radius must be positive");
    }
}

}  // namespace

NeighborSearch::NeighborSearch(double support_radius, Vector3 origin)
    : support_radius_(support_radius),
      origin_(origin),
      cells_(origin, support_radius) {
    validateSupportRadius(support_radius_);
}

TypedNeighborList NeighborSearch::buildTypedNeighborList(const ParticleSet& particles) {
    cells_.build(particles);

    TypedNeighborList neighbors;
    neighbors.fluid.resize(particles.size());
    neighbors.wall.resize(particles.size());
    const double support_radius_squared = support_radius_ * support_radius_;

    for (std::size_t index = 0; index < particles.size(); ++index) {
        const Vector3& position = particles.positions()[index];
        const CellIndex center_cell = cells_.cellIndex(position);

        for (int dz = -1; dz <= 1; ++dz) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dx = -1; dx <= 1; ++dx) {
                    const CellIndex candidate_cell{
                        center_cell.x + dx,
                        center_cell.y + dy,
                        center_cell.z + dz,
                    };

                    for (const std::size_t candidate : cells_.particlesInCell(candidate_cell)) {
                        if (candidate == index) {
                            continue;
                        }

                        const Vector3 offset = particles.positions()[candidate] - position;
                        if (normSquared(offset) <= support_radius_squared) {
                            if (particles.types()[candidate] == ParticleType::Fluid) {
                                neighbors.fluid[index].push_back(candidate);
                            } else if (particles.types()[candidate] == ParticleType::Wall) {
                                neighbors.wall[index].push_back(candidate);
                            }
                        }
                    }
                }
            }
        }

        std::sort(neighbors.fluid[index].begin(), neighbors.fluid[index].end());
        std::sort(neighbors.wall[index].begin(), neighbors.wall[index].end());
    }

    return neighbors;
}

NeighborList NeighborSearch::buildNeighborList(const ParticleSet& particles) {
    return combineNeighborList(buildTypedNeighborList(particles));
}

NeighborList NeighborSearch::combineNeighborList(const TypedNeighborList& neighbors) const {
    if (neighbors.fluid.size() != neighbors.wall.size()) {
        throw std::runtime_error("Typed neighbor list sizes must match");
    }

    NeighborList combined(neighbors.fluid.size());
    for (std::size_t index = 0; index < neighbors.fluid.size(); ++index) {
        combined[index].reserve(neighbors.fluid[index].size() + neighbors.wall[index].size());
        combined[index].insert(combined[index].end(), neighbors.fluid[index].begin(), neighbors.fluid[index].end());
        combined[index].insert(combined[index].end(), neighbors.wall[index].begin(), neighbors.wall[index].end());
        std::sort(combined[index].begin(), combined[index].end());
    }

    return combined;
}

std::vector<NeighborCountSummary> NeighborSearch::countNeighborsByType(const TypedNeighborList& neighbors) const {
    if (neighbors.fluid.size() != neighbors.wall.size()) {
        throw std::runtime_error("Typed neighbor list sizes must match");
    }

    std::vector<NeighborCountSummary> summaries(neighbors.fluid.size());
    for (std::size_t index = 0; index < neighbors.fluid.size(); ++index) {
        summaries[index].fluid = neighbors.fluid[index].size();
        summaries[index].wall = neighbors.wall[index].size();
        summaries[index].total = summaries[index].fluid + summaries[index].wall;
    }

    return summaries;
}

std::vector<NeighborCountSummary> NeighborSearch::countNeighborsByType(
    const ParticleSet& particles,
    const NeighborList& neighbors) const {
    if (neighbors.size() != particles.size()) {
        throw std::runtime_error("Neighbor list size must match particle count");
    }

    std::vector<NeighborCountSummary> summaries(neighbors.size());
    for (std::size_t index = 0; index < neighbors.size(); ++index) {
        summaries[index].total = neighbors[index].size();
        for (const std::size_t neighbor : neighbors[index]) {
            if (neighbor >= particles.size()) {
                throw std::runtime_error("Neighbor index out of range");
            }

            if (particles.types()[neighbor] == ParticleType::Fluid) {
                ++summaries[index].fluid;
            } else if (particles.types()[neighbor] == ParticleType::Wall) {
                ++summaries[index].wall;
            }
        }
    }

    return summaries;
}

void NeighborSearch::updateNeighborCounts(ParticleSet& particles, const TypedNeighborList& neighbors) const {
    if (neighbors.fluid.size() != particles.size() || neighbors.wall.size() != particles.size()) {
        throw std::runtime_error("Typed neighbor list size must match particle count");
    }

    const std::vector<NeighborCountSummary> summaries = countNeighborsByType(neighbors);

    std::vector<std::size_t>& total_counts = particles.neighborCounts();
    std::vector<std::size_t>& fluid_counts = particles.fluidNeighborCounts();
    std::vector<std::size_t>& wall_counts = particles.wallNeighborCounts();
    for (std::size_t index = 0; index < summaries.size(); ++index) {
        total_counts[index] = summaries[index].total;
        fluid_counts[index] = summaries[index].fluid;
        wall_counts[index] = summaries[index].wall;
    }
}

void NeighborSearch::updateNeighborCounts(ParticleSet& particles, const NeighborList& neighbors) const {
    const std::vector<NeighborCountSummary> summaries = countNeighborsByType(particles, neighbors);

    std::vector<std::size_t>& total_counts = particles.neighborCounts();
    std::vector<std::size_t>& fluid_counts = particles.fluidNeighborCounts();
    std::vector<std::size_t>& wall_counts = particles.wallNeighborCounts();
    for (std::size_t index = 0; index < summaries.size(); ++index) {
        total_counts[index] = summaries[index].total;
        fluid_counts[index] = summaries[index].fluid;
        wall_counts[index] = summaries[index].wall;
    }
}

double NeighborSearch::supportRadius() const noexcept {
    return support_radius_;
}

const Vector3& NeighborSearch::origin() const noexcept {
    return origin_;
}

const CellLinkedList& NeighborSearch::cells() const noexcept {
    return cells_;
}

}  // namespace lsmps
