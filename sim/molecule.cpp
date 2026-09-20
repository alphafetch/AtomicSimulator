#include "molecule.hpp"

#include "bondFind.hpp"
#include "../core/vector.hpp"

std::unordered_map<int, Molecule> buildMolecules(std::vector<Atom>& atoms, std::vector<int>& parent) {
    std::unordered_map<int, Molecule> molecules;
    for (auto& atom : atoms) {
        int root = find(parent, atom.id);

        molecules[root].mass += atom.mass;

        molecules[root].composition[atom.element]++;

        molecules[root].posSum.x += atom.pos.x;
        molecules[root].posSum.y += atom.pos.y;

        molecules[root].atomCount++;
    }

    for (auto& molecule : molecules) {
        auto& m = molecule.second;

        std::string compStr;
        for (auto& elem : m.composition) {
            compStr += elem.first;
            int count = elem.second;
            if (elem.second != 1) {
                compStr += std::to_string(count);
            }
        }

        m.formula = compStr;

        m.centeroid = vec::Vector2(m.posSum.x / m.atomCount, m.posSum.y / m.atomCount);
    }

    return molecules;
}