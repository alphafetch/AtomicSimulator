#ifndef SIM_PARTICLE_HPP
#define SIM_PARTICLE_HPP

#include <string>
#include <array>

#include "../core/vector.hpp"

struct Atom {
    int id;

    // Individual particle stats
    vec::Vector3 pos;
    vec::Vector3 vel;

    // Mass in kilograms
    double mass;
    
    // What particle is this?
    int electrons;
    int protons;
    int neutrons;

    // Element symbol (e.g. "He"), set via the periodic table lookup after
    // construction. "?" means no match was found (e.g. 0 protons).
    std::string element = "?";
    std::array<int, 3> subatomTier;

    Atom(float x, float y, float z, float vx, float vy, float vz, int electrons, int protons, int neutrons, int id);
    std::array<int, 3> getSubatomTier();
};

#endif