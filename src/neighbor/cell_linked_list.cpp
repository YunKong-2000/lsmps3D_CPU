#include "neighbor/cell_linked_list.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <stdexcept>

namespace lsmps {
namespace {

int coordinateToCell(double coordinate, double origin, double cell_size) {
    return static_cast<int>(std::floor((coordinate - origin) / cell_size));
}

void validateCellSize(double cell_size) {
    if (cell_size <= 0.0) {
        throw std::runtime_error("CellLinkedList cell size must be positive");
    }
}

}  // namespace

std::size_t CellIndexHash::operator()(const CellIndex& index) const noexcept {
    const std::size_t hx = std::hash<int>{}(index.x);
    const std::size_t hy = std::hash<int>{}(index.y);
    const std::size_t hz = std::hash<int>{}(index.z);
    return hx ^ (hy + 0x9e3779b9U + (hx << 6U) + (hx >> 2U)) ^
           (hz + 0x9e3779b9U + (hy << 6U) + (hy >> 2U));
}

CellLinkedList::CellLinkedList(Vector3 origin, double cell_size)
    : origin_(origin), cell_size_(cell_size) {
    validateCellSize(cell_size_);
}

void CellLinkedList::reset(Vector3 origin, double cell_size) {
    validateCellSize(cell_size);
    origin_ = origin;
    cell_size_ = cell_size;
    clear();
}

void CellLinkedList::build(const ParticleSet& particles) {
    clear();
    for (std::size_t index = 0; index < particles.size(); ++index) {
        cells_[cellIndex(particles.positions()[index])].push_back(index);
    }
}

void CellLinkedList::clear() {
    cells_.clear();
}

CellIndex CellLinkedList::cellIndex(const Vector3& position) const {
    return {
        coordinateToCell(position.x, origin_.x, cell_size_),
        coordinateToCell(position.y, origin_.y, cell_size_),
        coordinateToCell(position.z, origin_.z, cell_size_),
    };
}

const std::vector<std::size_t>& CellLinkedList::particlesInCell(const CellIndex& index) const {
    const auto iter = cells_.find(index);
    if (iter == cells_.end()) {
        return empty_cell_;
    }
    return iter->second;
}

const Vector3& CellLinkedList::origin() const noexcept {
    return origin_;
}

double CellLinkedList::cellSize() const noexcept {
    return cell_size_;
}

std::size_t CellLinkedList::cellCount() const noexcept {
    return cells_.size();
}

}  // namespace lsmps
