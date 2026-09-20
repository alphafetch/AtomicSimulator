#ifndef SIM_UTILS_HPP
#define SIM_UTILS_HPP

#include "atom.hpp"
#include "../core/vector.hpp"

vec::Vector2 getPointOfCenter(vec::Vector2 center, float radius, float rad);

float getAtomicDistance(Atom a, Atom b);

vec::Vector2 getForceMagnitude(float dist, Atom a, Atom b);

#endif