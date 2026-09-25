#ifndef SIM_UTILS_HPP
#define SIM_UTILS_HPP

#include <filesystem>

#include "atom.hpp"
#include "../core/vector.hpp"
#include "../core/settings.hpp"

vec::Vector2 getPointOfCenter(vec::Vector2 center, float radius, float rad);

float getAtomicDistance(Atom a, Atom b);

vec::Vector2 getForceMagnitude(float dist, Atom a, Atom b, const Settings& settings);

std::filesystem::path getUserHomeDir();

#endif