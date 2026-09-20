#include "vector.hpp"

vec::Vector1 vec::Vector1::operator*(const vec::Vector1& other) {
    Vector1 v(this->x * other.x);
    return v;
}

vec::Vector1 vec::Vector1::operator*(const float other) {
    Vector1 v(this->x * other);
    return v;
}

vec::Vector2 vec::Vector2::operator*(const vec::Vector2& other) {
    Vector2 v(this->x * other.x, this->y * other.y);
    return v;
}

vec::Vector2 vec::Vector2::operator*(const float other) {
    Vector2 v(this->x * other, this->y * other);
    return v;
}

vec::Vector3 vec::Vector3::operator*(const vec::Vector3& other) {
    Vector3 v(this->x * other.x, this->y * other.y, this->z * other.z);
    return v;
}

vec::Vector3 vec::Vector3::operator*(const float other) {
    Vector3 v(this->x * other, this->y * other, this->z * other);
    return v;
}