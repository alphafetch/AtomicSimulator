#ifndef SIM_PARTICLE_HPP
#define SIM_PARTICLE_HPP

#include <string>
#include <array>
#include <random>
#include <unordered_map>

#include "../core/vector.hpp"
#include "../core/elements.hpp"

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

    // Atom marked for removal?
    bool markedForRemoval = false;

    // Element symbol (e.g. "He"), set via the periodic table lookup after
    // construction. "?" means no match was found (e.g. 0 protons).
    std::string element = "?";
    std::array<int, 3> subatomTier;

    void refreshDerivedFields(const std::unordered_map<int, Element>& elementTable);
    Atom(float x, float y, float z, float vx, float vy, float vz, int electrons, int protons, int neutrons, int id, std::unordered_map<int, Element>& eT);
    std::array<int, 3> getSubatomTier();
};

void applyDecay(std::mt19937_64& rng, Atom& atom, std::unordered_map<int, Element>& elementTable);

#endif