#pragma once

#include "core/particle_set.hpp"
#include "core/simulation_config.hpp"
#include "neighbor/neighbor_search.hpp"

#include <cstddef>
#include <vector>

namespace lsmps {

struct FreeSurfaceDiagnostics {
    std::vector<double> open_ratio;
    std::vector<double> cone_ratio;
    std::vector<double> accessible_area_ratio;
    std::vector<int> reason_code;
};

class FreeSurfaceDetector {
public:
    explicit FreeSurfaceDetector(const FreeSurfaceConfig& config, double particle_spacing);

    FreeSurfaceDiagnostics detect(ParticleSet& particles, const TypedNeighborList& neighbors) const;

    const std::vector<Vector3>& directions() const noexcept;
    const std::vector<double>& directionWeights() const noexcept;

private:
    FreeSurfaceConfig config_;
    double particle_spacing_ = 1.0;
    double screen_radius_ = 1.0;
    double wall_patch_radius_ = 1.0;
    double particle_radius_ = 0.5;
    double cone_cosine_ = 0.0;

    std::vector<Vector3> directions_;
    std::vector<double> direction_weights_;
    std::vector<std::vector<std::size_t>> cone_direction_indices_;
    std::vector<double> cone_areas_;

    void buildCubedSphere(std::size_t q);
    void buildConeTable();
};

}  // namespace lsmps
