#include "atom.hpp"

#include <iostream>
#include <cstdlib>

#include "../core/constants.hpp"

Atom::Atom(float x, float y, float z, float vx, float vy, float vz, int electrons, int protons, int neutrons, int id) {
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
    this->subatomTier = this->getSubatomTier();

    this->mass = 
        (electrons * ELECTRON_MASS_KG) + 
        (protons * PROTON_MASS_KG) + 
        (neutrons * NEUTRON_MASS_KG);
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