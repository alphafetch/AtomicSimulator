#ifndef SIM_PROTEIN_BUILDER_HPP
#define SIM_PROTEIN_BUILDER_HPP

#include <vector>
#include <unordered_map>

#include "molecule.hpp"

struct ProteinBuilderStorage {
    int activeElement = 0;
    int selectedAtomIndex = -1;
    std::vector<std::pair<size_t, size_t>> bonds;
    std::unordered_map<int, Molecule> molecules;
    std::vector<int> parent;
    bool isSaving = false;
    char saveBuf[128];

    void resizeParent(size_t s);
};

#endif