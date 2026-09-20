#ifndef CORE_UTILS_HPP
#define CORE_UTILS_HPP

#include <random>

#include "vector.hpp"

vec::Vector2 randomAtomPos(std::mt19937_64& rng);

vec::Vector2 randomAtomVel(std::mt19937_64& rng);

float rad(float deg);

#endif