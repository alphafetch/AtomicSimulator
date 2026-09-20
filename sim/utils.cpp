#include "utils.hpp"

#include <cmath>
#include <algorithm>

#include "../core/constants.hpp"

vec::Vector2 getPointOfCenter(vec::Vector2 center, float radius, float rad) {
    float x = center.x + radius * cosf(rad);
    float y = center.y + radius * sinf(rad);

    vec::Vector2 v;
    v.x = x;
    v.y = y;

    return v;
}

float getAtomicDistance(Atom a, Atom b) {
    return 
        sqrt(((b.pos.x - a.pos.x) * (b.pos.x - a.pos.x)) 
        + ((b.pos.y - a.pos.y) * (b.pos.y - a.pos.y)));
}

vec::Vector2 getForceMagnitude(float dist, Atom a, Atom b) {
    if (dist == 0) {
        dist = 1.0f;
    }

    float dx = b.pos.x - a.pos.x;
    float dy = b.pos.y - a.pos.y;
    vec::Vector2 dir(dx / dist, dy / dist);

    if (dist < FORCE_MIN_DIST) {
        return dir * FORCE_STRENGTH;
    } else if (dist > FORCE_MIN_DIST && dist < FORCE_MAX_DIST) {
        return (dir * -1.0f) * FORCE_STRENGTH;
    } else {
        return vec::Vector2(0.0f, 0.0f);
    }
}