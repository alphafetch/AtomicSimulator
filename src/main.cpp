#include <random>
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib>
#include <numbers>
#include <vector>
#include <array>
#include <cmath>
#include <chrono>

#include "raylib.h"
#define RAYGUI_IMPLEMENTATION
#include "../external/raygui.h"

#include "../core/constants.hpp"
#include "../core/utils.hpp"
#include "../core/vector.hpp"
#include "../core/elements.hpp"

#include "../sim/atom.hpp"
#include "../sim/utils.hpp"
#include "../sim/bondFind.hpp"
#include "../sim/molecule.hpp"

#include "../core/settings.hpp"

using std::cout;

int main() {
    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, WINDOW_TITLE);
    SetTargetFPS(60);
    auto elementTable = loadElementTable("data/elements.csv");
    Settings settings;

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
        Atom atom(pos.x, pos.y, 0.0f, vel.x, vel.y, 0.0f, pct, pct, nct, i, elementTable);
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
    bool showSettings = false;
    std::vector<int> parent(atoms.size());
    std::vector<std::pair<size_t, size_t>> toFuse;
    std::vector<bool> fusing(atoms.size(), false);

    while (!WindowShouldClose()) {
        std::chrono::high_resolution_clock::time_point loopStartProfiler;
        if (DEBUG_PROFILER) loopStartProfiler = std::chrono::high_resolution_clock::now();

        // Reset the parent array
        for (size_t i = 0; i < parent.size(); i++) {
            parent[i] = i;
        }
        
        for (size_t i = 0; i < fusing.size(); i++) {
            fusing[i] = false;
        }

        // Change position based on atomic velocity
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i].pos.x += atoms[i].vel.x;
            atoms[i].pos.y += atoms[i].vel.y;
            atoms[i].pos.z += atoms[i].vel.z;

            if (atoms[i].pos.x >= WINDOW_WIDTH - (ELECTRON_FLOAT_RADIUS + ELECTRON_RENDER_RADIUS) || atoms[i].pos.x <= 0) {
                atoms[i].vel.x *= -1;
            }

            if (atoms[i].pos.y >= WINDOW_HEIGHT - (ELECTRON_FLOAT_RADIUS + ELECTRON_RENDER_RADIUS) || atoms[i].pos.y <= 0) {
                atoms[i].vel.y *= -1;
            }
        }

        std::chrono::high_resolution_clock::time_point startN2Profiler;
        if (DEBUG_PROFILER) startN2Profiler = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < atoms.size(); i++) {
            for (size_t j = i + 1; j < atoms.size(); j++) {
                float dist = getAtomicDistance(atoms[i], atoms[j]);

                float safeTemp = std::max(settings.temp, 0.1f);
                float effectiveBondDist = settings.bondDist / safeTemp;
                if (dist < effectiveBondDist) {
                    bond(parent, i, j);
                    float effectiveFusionMax = std::min(settings.fusionProtonMax * settings.temp, 118.0f);
                    if (atoms[i].protons + atoms[j].protons < effectiveFusionMax
                        && (fusing[i] == false && fusing[j] == false)) {
                        toFuse.push_back(std::pair<size_t, size_t>{i, j});
                        fusing[i] = true;
                        fusing[j] = true;
                    }
                }

                vec::Vector2 mag = getForceMagnitude(dist, atoms[i], atoms[j], settings);

                atoms[i].vel.x -= mag.x;
                atoms[i].vel.y -= mag.y;

                atoms[j].vel.x += mag.x;
                atoms[j].vel.y += mag.y;
            }
        }
        if (DEBUG_PROFILER) {
            auto endN2Profiler = std::chrono::high_resolution_clock::now();
            auto microsecondsN2Profiler = std::chrono::duration_cast<std::chrono::microseconds>(endN2Profiler - startN2Profiler).count();
            if (frame % 60 == 0) cout << "PROFILER: N2: " << microsecondsN2Profiler<< "\n";
        }

        for (size_t i = 0; i < toFuse.size(); i++) {
            atoms[toFuse[i].first].protons += atoms[toFuse[i].second].protons;
            atoms[toFuse[i].first].neutrons += atoms[toFuse[i].second].neutrons;
            atoms[toFuse[i].second].markedForRemoval = true;
        }

        for (int i = static_cast<int>(atoms.size() - 1); i >= 0; i--) {
            if (atoms[i].markedForRemoval) {
                atoms.erase(atoms.begin() + i);
            }
        }

        toFuse.clear();
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i].refreshDerivedFields(elementTable);
        }

        // Apply friction and decay to each atom
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i].vel = atoms[i].vel * settings.atomicFriction;
            applyDecay(rng, atoms[i], elementTable, settings);
        }
        
        auto molecules = buildMolecules(atoms, parent);
        
        // Log to a CSV and Console if needed
        for (size_t i = 0; i < atoms.size(); i++) {
            if (DEBUG_CONSOLE == 1 || DEBUG_CONSOLE == 3) { 
                cout 
                    << "Atom " << atoms[i].id << " - " << atoms[i].element
                    << ": X, Y, Z = " << atoms[i].pos.x << ", " << atoms[i].pos.y << ", " << atoms[i].pos.z << 
                    " - VX, VY, VZ = " << atoms[i].vel.x << ", " << atoms[i].vel.y << ", " << atoms[i].vel.z 
                    << " - Mass: " << atoms[i].mass
                    << "\n";
            }

            if (DEBUG_CONSOLE == 2 || DEBUG_CONSOLE == 3) {
                cout
                    << atoms[i].id << ": "
                    << find(parent, i) << "\n";
            }

            if (DEBUG_CSV) {
                Molecule& molecule = molecules[find(parent, i)];
                csvstream
                    << frame << ","
                    << atoms[i].id << ","
                    << atoms[i].pos.x << "," << atoms[i].pos.y << "," << atoms[i].pos.z << ","
                    << atoms[i].vel.x << "," << atoms[i].vel.y << "," << atoms[i].vel.z << ","
                    << atoms[i].mass << ","
                    << atoms[i].element << ","
                    << find(parent, i) << ","
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

        for (size_t i = 0; i < atoms.size(); i++) {
            DrawCircleLinesV(Vector2({atoms[i].pos.x, atoms[i].pos.y}), ELECTRON_FLOAT_RADIUS, Fade(RED, 0.5f));

            // Add labels to each element
            Molecule& molecule = molecules[find(parent, i)];
            if (molecule.atomCount == 1) {
                const char* sym = atoms[i].element.c_str();
                int textWidth = MeasureText(sym, (int)settings.fontSize);
                DrawText(sym, atoms[i].pos.x - (textWidth / 2), atoms[i].pos.y + RELATIVE_TEXT_HEIGHT_SOLO, settings.fontSize, Fade(GRAY, 0.8f));
            } else if (!molecule.labelDrawn) {
                std::string form = molecule.formula;
                std::string label = "[ " + form + "]";
                const char* moleculeLabel = label.c_str();
                int textWidth = MeasureText(moleculeLabel, (int)settings.fontSize);
                DrawText(moleculeLabel, molecule.centeroid.x - (textWidth / 2), molecule.centeroid.y + RELATIVE_TEXT_HEIGHT_MOLECULE, settings.fontSize, Fade(GRAY, 0.8f));
                molecule.labelDrawn = true;
            }
        
            for (int j = 0; j < atoms[i].subatomTier[2]; j++) {
                float theta = j * ((2 * std::numbers::pi) / atoms[i].subatomTier[2]) + frame * settings.electronOrbitSpeed;

                vec::Vector2 pt = getPointOfCenter(vec::Vector2(atoms[i].pos.x, atoms[i].pos.y), ELECTRON_FLOAT_RADIUS, theta);

                DrawCircleV(Vector2({pt.x, pt.y}), ELECTRON_RENDER_RADIUS, RED);
            }

            for (int j = 0; j < atoms[i].subatomTier[1] + atoms[i].subatomTier[0]; j++) {
                float theta = j * ((2 * std::numbers::pi) / (atoms[i].subatomTier[1] + atoms[i].subatomTier[0]));

                vec::Vector2 pt = getPointOfCenter(vec::Vector2(atoms[i].pos.x, atoms[i].pos.y), NUCLEAR_FLOAT_RADIUS, theta);

                if (j < atoms[i].subatomTier[1]) { DrawCircleV(Vector2({pt.x, pt.y}), SUBATOMIC_RENDER_RADIUS, GRAY); }
                else { DrawCircleV(Vector2({pt.x, pt.y}), SUBATOMIC_RENDER_RADIUS, GREEN); }
            }
        }

        if (showSettings) {
            int x = GuiWindowBox({0, 0, 600, 600}, "Settings");
            if (x) showSettings = !showSettings;
            GuiSlider(
                {150, 40, 300, 20}, 
                "Atomic Friction", 
                TextFormat("%.3f", settings.atomicFriction), 
                &settings.atomicFriction, 0.90f, 0.9999f
            );
            GuiSlider(
                {150, 70, 300, 20}, 
                "Force Strength", 
                TextFormat("%.3f", settings.forceStrength), 
                &settings.forceStrength, 0.0f, 1.0f
            );
            GuiSlider(
                {150, 100, 300, 20}, 
                "Force Minimum Distance", 
                TextFormat("%.f", settings.forceMinDist), 
                &settings.forceMinDist, 1.0f, 40.0f
            );
            GuiSlider(
                {150, 130, 300, 20}, 
                "Force Maximum Distance", 
                TextFormat("%.f", settings.forceMaxDist), 
                &settings.forceMaxDist, 10.0f, 200.0f
            );
            GuiSlider(
                {150, 160, 300, 20}, 
                "Bond Distance", 
                TextFormat("%.f", settings.bondDist), 
                &settings.bondDist, 5.0f, 40.0f
            );
            GuiSlider(
                {150, 190, 300, 20}, 
                "Decay Chance", 
                TextFormat("%.3f", settings.decayChance), 
                &settings.decayChance, 0.0f, 0.5f
            );
            GuiSlider(
                {150, 220, 300, 20}, 
                "Decay Proton Threshold", 
                TextFormat("%.f", settings.decayProtonThreshold), 
                &settings.decayProtonThreshold, 1.0f, 118.0f
            );
            GuiSlider(
                {150, 250, 300, 20}, 
                "Fusion Proton Max", 
                TextFormat("%.f", settings.fusionProtonMax), 
                &settings.fusionProtonMax, 2.0f, 118.0f
            );
            GuiSlider(
                {150, 280, 300, 20}, 
                "Electron Orbit Speed", 
                TextFormat("%.3f", settings.electronOrbitSpeed), 
                &settings.electronOrbitSpeed, 0.0f, 0.1f
            );
            GuiSlider(
                {150, 310, 300, 20}, 
                "Font Size", 
                TextFormat("%.f", settings.fontSize), 
                &settings.fontSize, 4.0f, 14.0f
            );
            GuiSlider(
                {150, 340, 300, 20}, 
                "Temperature", 
                TextFormat("%.3f", settings.temp), 
                &settings.temp, 0.0f, 3.0f
            );
        }

        if (GuiButton({710, 560, 70, 20}, "Settings")) {
            showSettings = !showSettings;
        }

        EndDrawing();

        if (DEBUG_PROFILER) {
            auto endLoopProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsLoopProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endLoopProfiler - loopStartProfiler).count();
            if (frame % 60 == 0) cout << "PROFILER: LOOP: " << microsecondsLoopProfiler << "\n";
        }

        frame++;
    }
    
    CloseWindow();
}
