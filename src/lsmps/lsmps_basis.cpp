#include "lsmps/lsmps_basis.hpp"

namespace lsmps {

LsmpsBasisVector evaluateTypeABasis(const Vector3& offset, double support_radius) {
    const double inv_re = 1.0 / support_radius;
    const double inv_re2 = inv_re * inv_re;
    const double dx = offset.x;
    const double dy = offset.y;
    const double dz = offset.z;

    LsmpsBasisVector basis;
    basis << dx * inv_re,
        dy * inv_re,
        dz * inv_re,
        0.5 * dx * dx * inv_re2,
        0.5 * dy * dy * inv_re2,
        0.5 * dz * dz * inv_re2,
        dx * dy * inv_re2,
        dx * dz * inv_re2,
        dy * dz * inv_re2;
    return basis;
}

LsmpsBasisVector evaluateTypeANeumannBasis(
    const Vector3& offset,
    const Vector3& wall_normal,
    double support_radius) {
    const double inv_re = 1.0 / support_radius;
    const double dx = offset.x;
    const double dy = offset.y;
    const double dz = offset.z;
    const double nx = wall_normal.x;
    const double ny = wall_normal.y;
    const double nz = wall_normal.z;

    LsmpsBasisVector basis;
    basis << nx,
        ny,
        nz,
        nx * dx * inv_re,
        ny * dy * inv_re,
        nz * dz * inv_re,
        (nx * dy + ny * dx) * inv_re,
        (nx * dz + nz * dx) * inv_re,
        (ny * dz + nz * dy) * inv_re;
    return basis;
}

}  // namespace lsmps
