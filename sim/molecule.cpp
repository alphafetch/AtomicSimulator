#include "molecule.hpp"

#include <set>
#include <string>

#include "bondFind.hpp"
#include "../core/vector.hpp"
#include "../core/constants.hpp"

bool decideProtein(const std::map<std::string, int>& comp, int atomCount) {
    std::set<std::string> chnops = {"C", "H", "N", "O", "P", "S"};
    if (atomCount >= MIN_PROTEIN_SIZE) {
        for (auto atom : comp) {
            if (chnops.count(atom.first)) {
                continue;
            } else {
                return false;
            }
        }

        return true;
    }

    return false;
}

std::unordered_map<int, Molecule> buildMolecules(std::vector<Atom>& atoms, std::vector<int>& parent) {
    std::unordered_map<int, Molecule> molecules;
    for (size_t i = 0; i < atoms.size(); i++) {
        int root = find(parent, i);

        molecules[root].mass += atoms[i].mass;
        molecules[root].maxY = atoms[i].pos.y < molecules[root].maxY ? atoms[i].pos.y : molecules[root].maxY;

        molecules[root].composition[atoms[i].element]++;

        molecules[root].posSum.x += atoms[i].pos.x;
        molecules[root].posSum.y += atoms[i].pos.y;

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
            compStr += " ";
        }

        m.formula = compStr;
        m.isProtein = decideProtein(m.composition, m.atomCount);

        m.centeroid = vec::Vector2(m.posSum.x / m.atomCount, m.posSum.y / m.atomCount);
    }

    return molecules;
}