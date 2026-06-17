#pragma once

#include <cmath>

namespace lsmps {

struct Vector3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    constexpr Vector3() = default;
    constexpr Vector3(double x_value, double y_value, double z_value)
        : x(x_value), y(y_value), z(z_value) {}

    constexpr Vector3& operator+=(const Vector3& rhs) {
        x += rhs.x;
        y += rhs.y;
        z += rhs.z;
        return *this;
    }

    constexpr Vector3& operator-=(const Vector3& rhs) {
        x -= rhs.x;
        y -= rhs.y;
        z -= rhs.z;
        return *this;
    }

    constexpr Vector3& operator*=(double scalar) {
        x *= scalar;
        y *= scalar;
        z *= scalar;
        return *this;
    }

    constexpr Vector3& operator/=(double scalar) {
        x /= scalar;
        y /= scalar;
        z /= scalar;
        return *this;
    }
};

constexpr Vector3 operator+(Vector3 lhs, const Vector3& rhs) {
    lhs += rhs;
    return lhs;
}

constexpr Vector3 operator-(Vector3 lhs, const Vector3& rhs) {
    lhs -= rhs;
    return lhs;
}

constexpr Vector3 operator-(const Vector3& value) {
    return {-value.x, -value.y, -value.z};
}

constexpr Vector3 operator*(Vector3 value, double scalar) {
    value *= scalar;
    return value;
}

constexpr Vector3 operator*(double scalar, Vector3 value) {
    value *= scalar;
    return value;
}

constexpr Vector3 operator/(Vector3 value, double scalar) {
    value /= scalar;
    return value;
}

constexpr double dot(const Vector3& lhs, const Vector3& rhs) {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

constexpr Vector3 cross(const Vector3& lhs, const Vector3& rhs) {
    return {
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x,
    };
}

inline double normSquared(const Vector3& value) {
    return dot(value, value);
}

inline double norm(const Vector3& value) {
    return std::sqrt(normSquared(value));
}

}  // namespace lsmps
