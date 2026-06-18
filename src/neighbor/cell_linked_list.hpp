#pragma once

#include "core/particle_set.hpp"
#include "core/vector3.hpp"

#include <cstddef>
#include <unordered_map>
#include <vector>

namespace lsmps {

struct CellIndex {
    int x = 0;
    int y = 0;
    int z = 0;

    bool operator==(const CellIndex& rhs) const noexcept {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }
};

struct CellIndexHash {
    std::size_t operator()(const CellIndex& index) const noexcept;
};

class CellLinkedList {
public:
    CellLinkedList() = default;
    CellLinkedList(Vector3 origin, double cell_size);

    void reset(Vector3 origin, double cell_size);
    void build(const ParticleSet& particles);
    void clear();

    CellIndex cellIndex(const Vector3& position) const;
    const std::vector<std::size_t>& particlesInCell(const CellIndex& index) const;

    const Vector3& origin() const noexcept;
    double cellSize() const noexcept;
    std::size_t cellCount() const noexcept;

private:
    Vector3 origin_;
    double cell_size_ = 1.0;
    std::unordered_map<CellIndex, std::vector<std::size_t>, CellIndexHash> cells_;
    std::vector<std::size_t> empty_cell_;
};

}  // namespace lsmps
