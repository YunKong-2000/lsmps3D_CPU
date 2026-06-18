#pragma once

#include "core/particle_set.hpp"
#include "core/vector3.hpp"
#include "neighbor/cell_linked_list.hpp"

#include <cstddef>
#include <vector>

namespace lsmps {

using NeighborList = std::vector<std::vector<std::size_t>>;

struct TypedNeighborList {
    NeighborList fluid;
    NeighborList wall;
};

struct NeighborCountSummary {
    std::size_t total = 0;
    std::size_t fluid = 0;
    std::size_t wall = 0;
};

class NeighborSearch {
public:
    NeighborSearch(double support_radius, Vector3 origin = {});

    TypedNeighborList buildTypedNeighborList(const ParticleSet& particles);
    NeighborList buildNeighborList(const ParticleSet& particles);
    NeighborList combineNeighborList(const TypedNeighborList& neighbors) const;
    std::vector<NeighborCountSummary> countNeighborsByType(const TypedNeighborList& neighbors) const;
    std::vector<NeighborCountSummary> countNeighborsByType(
        const ParticleSet& particles,
        const NeighborList& neighbors) const;
    void updateNeighborCounts(ParticleSet& particles, const TypedNeighborList& neighbors) const;
    void updateNeighborCounts(ParticleSet& particles, const NeighborList& neighbors) const;

    double supportRadius() const noexcept;
    const Vector3& origin() const noexcept;
    const CellLinkedList& cells() const noexcept;

private:
    double support_radius_ = 1.0;
    Vector3 origin_;
    CellLinkedList cells_;
};

}  // namespace lsmps
