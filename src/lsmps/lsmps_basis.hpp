#pragma once

#include "core/vector3.hpp"

#include <Eigen/Dense>

namespace lsmps {

constexpr int lsmps_basis_size = 9;

using LsmpsBasisVector = Eigen::Matrix<double, lsmps_basis_size, 1>;
using LsmpsMomentMatrix = Eigen::Matrix<double, lsmps_basis_size, lsmps_basis_size>;

LsmpsBasisVector evaluateTypeABasis(const Vector3& offset, double support_radius);

LsmpsBasisVector evaluateTypeANeumannBasis(
    const Vector3& offset,
    const Vector3& wall_normal,
    double support_radius);

}  // namespace lsmps
