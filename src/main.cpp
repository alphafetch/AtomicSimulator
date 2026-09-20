#include <random>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <numbers>
#include <vector>
#include <array>
#include <cmath>

#include "raylib.h"

#include "../core/constants.hpp"
#include "../core/utils.hpp"
#include "../core/vector.hpp"
#include "../core/elements.hpp"

#include "../sim/atom.hpp"
#include "../sim/utils.hpp"
#include "../sim/bondFind.hpp"
#include "../sim/molecule.hpp"

using std::cout;

int main() {
    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, WINDOW_TITLE);
    SetTargetFPS(60);
    auto elementTable = loadElementTable("data/elements.csv");

    std::random_device rd;
    unsigned int seed = rd();
    std::mt19937_64 rng(seed);

    std::vector<double> weights(118);
    for (size_t i = 0; i < weights.size(); i++) {
        weights[i] = sqrt(static_cast<double>(elementTable.at(i + 1).abundancePpm));
    }

    std::discrete_distribution<int> proton_ct(weights.begin(), weights.end());

    std::vector<Atom> atoms; 

    // Initalize particles
    for (int i = 0; i < COUNT; i++) {
        int pct = proton_ct(rng) + 1;
        auto elem = elementTable.at(pct);
        int nct = elem.neutrons;

        vec::Vector2 pos = randomAtomPos(rng);
        vec::Vector2 vel = randomAtomVel(rng);
        Atom atom(pos.x, pos.y, 0.0f, vel.x, vel.y, 0.0f, pct, pct, nct, i);
        atom.element = elem.symbol;
        atoms.push_back(atom);
    }

    // Initialize logging
    std::ofstream csvstream;
    if (DEBUG_CSV) {
        csvstream.open("atoms-" + std::to_string(seed) + ".csv");
        if (!csvstream.is_open()) {
            std::cerr << "Failed to initialize logging at atoms-#.csv.\n";
            exit(EXIT_FAILURE);
        }

        csvstream << "# SEED: " << seed << "\n";
        csvstream << "frame,id,x,y,z,vx,vy,vz,mass,element,parent,molecule_formula,molecule_mass\n"; 
    }

    if (DEBUG_CONSOLE == 1 || DEBUG_CONSOLE == 3) { cout << "SEED: " << seed << "\n"; }

    // Other miscellaneous initializations
    int frame = 0;
    std::vector<int> parent(atoms.size());

    while (!WindowShouldClose()) {
        // Reset the parent array
        for (size_t i = 0; i < parent.size(); i++) {
            parent[i] = i;
        }

        // Change position based on atomic velocity
        for (auto& atom : atoms) {
            atom.pos.x += atom.vel.x;
            atom.pos.y += atom.vel.y;
            atom.pos.z += atom.vel.z;

            if (atom.pos.x >= WINDOW_WIDTH - (ELECTRON_FLOAT_RADIUS + ELECTRON_RENDER_RADIUS) || atom.pos.x <= 0) {
                atom.vel.x *= -1;
            }

            if (atom.pos.y >= WINDOW_HEIGHT - (ELECTRON_FLOAT_RADIUS + ELECTRON_RENDER_RADIUS) || atom.pos.y <= 0) {
                atom.vel.y *= -1;
            }
        }

        // Apply bounces
        for (size_t i = 0; i < atoms.size(); i++) {
            for (size_t j = i + 1; j < atoms.size(); j++) {
                float dist = getAtomicDistance(atoms[i], atoms[j]);

                if (dist < BOND_DIST) {
                    bond(parent, i, j);
                }

                vec::Vector2 mag = getForceMagnitude(dist, atoms[i], atoms[j]);

                atoms[i].vel.x -= mag.x;
                atoms[i].vel.y -= mag.y;

                atoms[j].vel.x += mag.x;
                atoms[j].vel.y += mag.y;
            }
        }
        
        auto molecules = buildMolecules(atoms, parent);

        // Apply friction to each atom
        for (auto& atom : atoms) {
            atom.vel = atom.vel * ATOMIC_FRICTION;
        }
        
        // Log to a CSV and Console if needed
        for (auto& atom : atoms) {
            if (DEBUG_CONSOLE == 1 || DEBUG_CONSOLE == 3) { 
                cout 
                    << "Atom " << atom.id << " - " << atom.element
                    << ": X, Y, Z = " << atom.pos.x << ", " << atom.pos.y << ", " << atom.pos.z << 
                    " - VX, VY, VZ = " << atom.vel.x << ", " << atom.vel.y << ", " << atom.vel.z 
                    << " - Mass: " << atom.mass
                    << "\n";
            }

            if (DEBUG_CONSOLE == 2 || DEBUG_CONSOLE == 3) {
                cout
                    << atom.id << ": "
                    << find(parent, atom.id) << "\n";
            }

            if (DEBUG_CSV) {
                Molecule& molecule = molecules[find(parent, atom.id)];
                csvstream
                    << frame << ","
                    << atom.id << ","
                    << atom.pos.x << "," << atom.pos.y << "," << atom.pos.z << ","
                    << atom.vel.x << "," << atom.vel.y << "," << atom.vel.z << ","
                    << atom.mass << ","
                    << atom.element << ","
                    << find(parent, atom.id) << ","
                    << molecule.formula << ","
                    << molecule.mass
                    << "\n";
            }
        }

        if (DEBUG_CONSOLE == 4 || DEBUG_CONSOLE == 3) {
            for (auto& entry : molecules) {
                Molecule& molecule = entry.second;
                cout
                    << "Molecule (root " << entry.first << "): " << molecule.formula
                    << " - Mass: " << molecule.mass
                    << "\n";
            }
        }

        // Draw to the screen (TODO)
        BeginDrawing();
        ClearBackground(BLACK);

        for (auto& atom : atoms) {
            DrawCircleLinesV(Vector2({atom.pos.x, atom.pos.y}), ELECTRON_FLOAT_RADIUS, Fade(RED, 0.5f));

            // Add labels to each element
            Molecule& molecule = molecules[find(parent, atom.id)];
            if (molecule.atomCount == 1) {
                const char* sym = atom.element.c_str();
                int textWidth = MeasureText(sym, FONT_SIZE);
                DrawText(sym, atom.pos.x - (textWidth / 2), atom.pos.y + RELATIVE_TEXT_HEIGHT_SOLO, FONT_SIZE, Fade(GRAY, 0.8f));
            } else if (!molecule.labelDrawn) {
                std::string form = molecule.formula;
                std::string label = "[ " + form + "]";
                const char* moleculeLabel = label.c_str();
                int textWidth = MeasureText(moleculeLabel, FONT_SIZE);
                DrawText(moleculeLabel, molecule.centeroid.x - (textWidth / 2), molecule.centeroid.y + RELATIVE_TEXT_HEIGHT_MOLECULE, FONT_SIZE, Fade(GRAY, 0.8f));
                molecule.labelDrawn = true;
            }
        
            for (int j = 0; j < atom.subatomTier[2]; j++) {
                float theta = j * ((2 * std::numbers::pi) / atom.subatomTier[2]) + frame * ELECTRON_ORBIT_SPEED;

                vec::Vector2 pt = getPointOfCenter(vec::Vector2(atom.pos.x, atom.pos.y), ELECTRON_FLOAT_RADIUS, theta);

                DrawCircleV(Vector2({pt.x, pt.y}), ELECTRON_RENDER_RADIUS, RED);
            }

            for (int j = 0; j < atom.subatomTier[1] + atom.subatomTier[0]; j++) {
                float theta = j * ((2 * std::numbers::pi) / (atom.subatomTier[1] + atom.subatomTier[0]));

                vec::Vector2 pt = getPointOfCenter(vec::Vector2(atom.pos.x, atom.pos.y), NUCLEAR_FLOAT_RADIUS, theta);

                if (j < atom.subatomTier[1]) { DrawCircleV(Vector2({pt.x, pt.y}), SUBATOMIC_RENDER_RADIUS, GRAY); }
                else { DrawCircleV(Vector2({pt.x, pt.y}), SUBATOMIC_RENDER_RADIUS, GREEN); }
            }
        }

        EndDrawing();

        frame++;
    }
    CloseWindow();
}
