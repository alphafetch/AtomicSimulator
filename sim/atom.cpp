#include "atom.hpp"

#include <iostream>
#include <cstdlib>

#include "../core/constants.hpp"

void Atom::refreshDerivedFields(const std::unordered_map<int, Element>& elementTable) { 
    this->subatomTier = this->getSubatomTier(); 
    
    this->electrons = this->protons;
    this->mass = (this->electrons * ELECTRON_MASS_KG) + 
                 (this->protons * PROTON_MASS_KG) + 
                 (this->neutrons * NEUTRON_MASS_KG); 

    auto elem = elementTable.at(this->protons);
    this->element = elem.symbol;
}

Atom::Atom(float x, float y, float z, float vx, float vy, float vz, int electrons, int protons, int neutrons, int id, std::unordered_map<int, Element>& eT) {
    this->id = id;

    this->pos.x = x;
    this->pos.y = y;
    this->pos.z = z;
    
    this->vel.x = vx;
    this->vel.y = vy;
    this->vel.z = vz;
    
    this->electrons = electrons;
    this->protons = protons;
    this->neutrons = neutrons;
    
    this->refreshDerivedFields(eT);
}

std::array<int, 3> Atom::getSubatomTier() {
    if (this->protons > 118 || this->protons < 1) {
        std::cerr << "Failed to get subatom tier.\n";
        exit(EXIT_FAILURE);
    }

    std::array<int, 3> arr;

    if (this->protons < 26) {
        arr[0] = 2;
        arr[1] = 1;
        arr[2] = 2;
    } else if (this->protons < 71) {
        arr[0] = 3;
        arr[1] = 3;
        arr[2] = 2;
    } else if (this->protons < 96) {
        arr[0] = 4;
        arr[1] = 4;
        arr[2] = 3;
    } else if (this->protons < 119) {
        arr[0] = 5;
        arr[1] = 5;
        arr[2] = 4;
    }

    return arr;
}

void applyDecay(std::mt19937_64& rng, Atom& atom, std::unordered_map<int, Element>& elementTable, const Settings& settings) {
    std::uniform_real_distribution<float> dist(0, 1);

    float roll = dist(rng);
    float effectiveDecayChance = settings.decayChance * settings.temp;
    if (roll <= effectiveDecayChance && atom.protons > settings.decayProtonThreshold) {
        atom.protons -= 2;
        atom.neutrons -= 2;
        atom.refreshDerivedFields(elementTable);
    }
}