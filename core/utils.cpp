#include "utils.hpp"

#include <numbers>

#include "constants.hpp"

using std::mt19937_64,
      std::uniform_real_distribution;

// Temporarily only generates X and Y, not yet Z
vec::Vector2 randomAtomPos(mt19937_64& rng) {
    uniform_real_distribution<float> dx(0, WINDOW_WIDTH);
    uniform_real_distribution<float> dy(0, WINDOW_HEIGHT);

    vec::Vector2 v;
    v.x = dx(rng);
    v.y = dy(rng);

    return v;
}

vec::Vector2 randomAtomVel(mt19937_64& rng) {
    uniform_real_distribution<float> dx(-1, 1);
    uniform_real_distribution<float> dy(-1, 1);

    vec::Vector2 v;
    v.x = dx(rng);
    v.y = dy(rng);

    return v;
}

float rad(float deg) { return deg * (std::numbers::pi / 180); }