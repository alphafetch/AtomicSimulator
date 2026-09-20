#ifndef SIM_MOLECULE_HPP
#define SIM_MOLECULE_HPP

#include <map>
#include <unordered_map>
#include <string>
#include <vector>

#include "atom.hpp"

struct Molecule {
    double mass = 0.0;
    std::map<std::string, int> composition;
    std::string formula;
    vec::Vector2 posSum;
    int atomCount = 0;
    vec::Vector2 centeroid;
    bool labelDrawn = false;
};

std::unordered_map<int, Molecule> buildMolecules(std::vector<Atom>& atoms, std::vector<int>& parent);
 
#endif