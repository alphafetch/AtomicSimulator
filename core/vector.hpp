#ifndef CORE_VECTOR_HPP
#define CORE_VECTOR_HPP

namespace vec {
    struct Vector1 {
        float x;

        Vector1() : x(0) {}
        Vector1(float x) : x(x) {}
        Vector1(const Vector1& other) = default;
        Vector1& operator=(const Vector1& other) = default;
        Vector1 operator*(const Vector1& other);
        Vector1 operator*(const float other);
    };

    struct Vector2 {
        float x;
        float y;

        Vector2() : x(0), y(0) {}
        Vector2(float x, float y) : x(x), y(y) {}
        Vector2(const Vector2& other) = default;
        Vector2& operator=(const Vector2& other) = default;
        Vector2 operator*(const Vector2& other);
        Vector2 operator*(const float other);
    };

    struct Vector3 {
        float x;
        float y;
        float z;
        
        Vector3() : x(0), y(0), z(0) {}
        Vector3(float x, float y, float z) : x(x), y(y), z(z) {}
        Vector3(const Vector3& other) = default;
        Vector3& operator=(const Vector3& other) = default;
        Vector3 operator*(const Vector3& other);
        Vector3 operator*(const float other);
    };
}

#endif