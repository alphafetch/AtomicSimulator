#include "screens.hpp"

#include <random>
#include <fstream>
#include <iostream>
#include <chrono>
#include <filesystem>

#include "raylib.h"
#include "../external/raygui.h"

#include "../core/constants.hpp"
#include "../core/elements.hpp"
#include "../core/settings.hpp"
#include "../core/utils.hpp"

#include "../sim/atom.hpp"
#include "../sim/utils.hpp"
#include "../sim/bondFind.hpp"
#include "../sim/molecule.hpp"
#include "../sim/render.hpp"
#include "../sim/proteinBuilder.hpp"

SimScreen runMainMenu() {
    SetWindowTitle(MAIN_WINDOW_TITLE);
    
    while (!WindowShouldClose()) {
        BeginDrawing();

        ClearBackground(DARKGRAY);

        if (GuiButton({340, 215, 120, 40}, "Blank Simulation")) {
            return BLANK_SIMULATION;
        }

        if (GuiButton({340, 265, 120, 40}, "Atom Simulation")) {
            return ATOM_SIMULATION;
        }

        if (GuiButton({340, 315, 120, 40}, "Protein Builder")) {
            return PROTEIN_BUILDER;
        }

        if (GuiButton({340, 365, 120, 40}, "Quit")) {
            break;
        }

        EndDrawing();
    }

    return QUIT;
}

SimScreen runAtomSimulator() {
    SetWindowTitle(ATOM_SIMULATOR_WINDOW_TITLE);

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

    if (DEBUG_CONSOLE == 1 || DEBUG_CONSOLE == 3) { std::cout << "SEED: " << seed << "\n"; }

    // Other miscellaneous initializations
    int frame = 0;
    bool showSettings = false;
    std::vector<int> parent(atoms.size());
    std::vector<std::pair<size_t, size_t>> toFuse;
    std::vector<bool> fusing(atoms.size(), false);

    while (!WindowShouldClose()) {
        std::chrono::high_resolution_clock::time_point loopStartProfiler;
        if (DEBUG_PROFILER) loopStartProfiler = std::chrono::high_resolution_clock::now();

        std::chrono::high_resolution_clock::time_point startResetLoopProfiler;
        if (DEBUG_PROFILER) startResetLoopProfiler = std::chrono::high_resolution_clock::now();
        // Reset the parent array
        for (size_t i = 0; i < parent.size(); i++) {
            parent[i] = i;
        }
        
        for (size_t i = 0; i < fusing.size(); i++) {
            fusing[i] = false;
        }
        if (DEBUG_PROFILER) {
            auto endResetLoopProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsResetLoopProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endResetLoopProfiler - startResetLoopProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: RESET LOOPS: " << microsecondsResetLoopProfiler << "\n";
        }

        std::chrono::high_resolution_clock::time_point startPosProfiler;
        if (DEBUG_PROFILER) startPosProfiler = std::chrono::high_resolution_clock::now();
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
        if (DEBUG_PROFILER) {
            auto endPosProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsPosProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endPosProfiler - startPosProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: POS: " << microsecondsPosProfiler << "\n";
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
            if (frame % 60 == 0) std::cout << "PROFILER: N2: " << microsecondsN2Profiler << "\n";
        }

        std::chrono::high_resolution_clock::time_point startFuseAppProfiler;
        if (DEBUG_PROFILER) startFuseAppProfiler = std::chrono::high_resolution_clock::now();
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

        if (DEBUG_PROFILER) {
            auto endFuseAppProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsFuseAppProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endFuseAppProfiler - startFuseAppProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: FUSE APP.: " << microsecondsFuseAppProfiler << "\n";
        }

        std::chrono::high_resolution_clock::time_point startRDFProfiler;
        if (DEBUG_PROFILER) startRDFProfiler = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i].refreshDerivedFields(elementTable);
        }
        if (DEBUG_PROFILER) {
            auto endRDFProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsRDFProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endRDFProfiler - startRDFProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: RDF: " << microsecondsRDFProfiler << "\n";
        }

        std::chrono::high_resolution_clock::time_point startFrictionDecayProfiler;
        if (DEBUG_PROFILER) startFrictionDecayProfiler = std::chrono::high_resolution_clock::now();
        // Apply friction and decay to each atom
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i].vel = atoms[i].vel * settings.atomicFriction;
            applyDecay(rng, atoms[i], elementTable, settings);
        }
        if (DEBUG_PROFILER) {
            auto endFrictionDecayProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsFrictionDecayProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endFrictionDecayProfiler - startFrictionDecayProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: FRIC/DECAY: " << microsecondsFrictionDecayProfiler << "\n";
        }
        
        std::chrono::high_resolution_clock::time_point startBuildMolecProfiler;
        if (DEBUG_PROFILER) startBuildMolecProfiler = std::chrono::high_resolution_clock::now();
        auto molecules = buildMolecules(atoms, parent);
        if (DEBUG_PROFILER) {
            auto endBuildMolecProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsBuildMolecProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endBuildMolecProfiler - startBuildMolecProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: BUILD MOLEC.: " << microsecondsBuildMolecProfiler << "\n";
        }
        
        // Log to a CSV and Console if needed
        for (size_t i = 0; i < atoms.size(); i++) {
            if (DEBUG_CONSOLE == 1 || DEBUG_CONSOLE == 3) { 
                std::cout 
                    << "Atom " << atoms[i].id << " - " << atoms[i].element
                    << ": X, Y, Z = " << atoms[i].pos.x << ", " << atoms[i].pos.y << ", " << atoms[i].pos.z << 
                    " - VX, VY, VZ = " << atoms[i].vel.x << ", " << atoms[i].vel.y << ", " << atoms[i].vel.z 
                    << " - Mass: " << atoms[i].mass
                    << "\n";
            }

            if (DEBUG_CONSOLE == 2 || DEBUG_CONSOLE == 3) {
                std::cout
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
                std::cout
                    << "Molecule (root " << entry.first << "): " << molecule.formula
                    << " - Mass: " << molecule.mass
                    << "\n";
            }
        }

        std::chrono::high_resolution_clock::time_point startDrawProfiler;
        if (DEBUG_PROFILER) startDrawProfiler = std::chrono::high_resolution_clock::now();
        // Draw to the screen (TODO)
        BeginDrawing();
        ClearBackground(BLACK);
        
        for (size_t i = 0; i < atoms.size(); i++) {
            renderAtom(atoms[i], settings, frame);
        }

        for (size_t i = 0; i < atoms.size(); i++) {
            // Add labels to each element
            Molecule& molecule = molecules[find(parent, i)];
            if (molecule.atomCount == 1) {
                const char* sym = atoms[i].element.c_str();
                drawCenteredLabel(atoms[i].pos.x, atoms[i].pos.y, sym, settings, RELATIVE_TEXT_HEIGHT_SOLO, Fade(GRAY, 0.8f));
            } else if (!molecule.labelDrawn) {
                std::string form = molecule.formula;
                std::string label = "[ " + form + "]";
                if (molecule.isProtein) label += " - Prtn.";
                const char* moleculeLabel = label.c_str();
                molecule.labelDrawn = true;
                drawCenteredLabel(molecule.centeroid.x, molecule.centeroid.y, moleculeLabel, settings, RELATIVE_TEXT_HEIGHT_MOLECULE, Fade(GRAY, 0.8f));
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
            if (settings.simpleAtoms) GuiDisable();
            GuiSlider(
                {150, 280, 300, 20}, 
                "Electron Orbit Speed", 
                TextFormat("%.3f", settings.electronOrbitSpeed), 
                &settings.electronOrbitSpeed, 0.0f, 0.1f
            );
            if (settings.simpleAtoms) GuiEnable();
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
            GuiCheckBox(
                {150, 370, 20, 20},
                "Simple Atoms (performance boost)",
                &settings.simpleAtoms
            );
        }

        if (GuiButton({710, 560, 70, 20}, "Settings")) {
            showSettings = !showSettings;
        }
        
        bool menu = GuiButton({10, 10, 30, 30}, "<");

        EndDrawing();

        if (menu) return MAIN;

        if (DEBUG_PROFILER) {
            auto endDrawProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsDrawProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endDrawProfiler - startDrawProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: DRAW: " << microsecondsDrawProfiler << "\n";
        }

        if (DEBUG_PROFILER) {
            auto endLoopProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsLoopProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endLoopProfiler - loopStartProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: LOOP: " << microsecondsLoopProfiler << "\n";
        }

        frame++;
    }
    
    return QUIT;
}

SimScreen runProteinBuilder() {
    SetWindowTitle(PROTEIN_BUILDER_WINDOW_TITLE);
    auto elementTable = loadElementTable("data/elements.csv");
    Settings settings;
    
    int frame = 0;
    std::vector<Atom> atoms;
    std::vector<int> protCounts = {1, 6, 7, 8, 15, 16};
    std::vector<int> neutCounts = {0, 6, 7, 8, 16, 16};
    ProteinBuilderStorage PBS;
    std::vector<Rectangle> UIBounds = {{0, 560, 500, 40}, {10, 10, 30, 30}, {700, 10, 90, 60}};
    bool isCreatingBond = false;

    while (!WindowShouldClose()) {

        for (size_t i = 0; i < atoms.size(); i++) {
            if (CheckCollisionPointCircle(GetMousePosition(), {atoms[i].pos.x, atoms[i].pos.y}, ELECTRON_FLOAT_RADIUS + ELECTRON_RENDER_RADIUS)
                && IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
                if (PBS.selectedAtomIndex == -1) PBS.selectedAtomIndex = i;
                else if (PBS.selectedAtomIndex == (int)i) {
                    PBS.selectedAtomIndex = -1;
                    isCreatingBond = true;
                } else {
                    PBS.bonds.push_back(std::pair<size_t, size_t>{PBS.selectedAtomIndex, i});
                    PBS.selectedAtomIndex = -1;
                    isCreatingBond = true;
                }
            }
        }
        
        for (size_t i = 0; i < PBS.parent.size(); i++) {
            PBS.parent[i] = i;
        }

        for (auto b : PBS.bonds) {
            bond(PBS.parent, b.first, b.second);
        }
        
        PBS.molecules = buildMolecules(atoms, PBS.parent);

        bool isUIClick = false;
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            for (auto& bound : UIBounds) {
                if (CheckCollisionPointRec(GetMousePosition(), bound)) {
                    isUIClick = true;
                    break;
                }
            }
        }

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) 
            && !isUIClick
            && PBS.selectedAtomIndex == -1
            && !isCreatingBond) {
            Vector2 pos = GetMousePosition();

            Atom atom(
                pos.x, pos.y, 0.0f,
                0.0f, 0.0f, 0.0f,
                protCounts[PBS.activeElement], 
                protCounts[PBS.activeElement], 
                neutCounts[PBS.activeElement],
                (int)atoms.size(),
                elementTable
            );

            atoms.push_back(atom);
            PBS.resizeParent(atoms.size());
        }

        if (isCreatingBond) isCreatingBond = false;

        BeginDrawing();

        ClearBackground(BLACK);

        for (auto bond : PBS.bonds) {
            DrawLineEx(
                {atoms[bond.first].pos.x, atoms[bond.first].pos.y}, 
                {atoms[bond.second].pos.x, atoms[bond.second].pos.y},
                BOND_THICKNESS,
                Fade(GRAY, 0.6f)
            );
        }

        for (size_t i = 0; i < atoms.size(); i++) {
            renderAtom(atoms[i], settings, frame);
            if (PBS.selectedAtomIndex == (int)i) {
                DrawCircleLinesV(
                    {atoms[i].pos.x, atoms[i].pos.y}, 
                    ELECTRON_FLOAT_RADIUS + ELECTRON_RENDER_RADIUS + SELECT_RING_ADDITIONAL_RAD, 
                    WHITE
                );
            }
        }

        for (size_t i = 0; i < atoms.size(); i++) {
            Molecule& molecule = PBS.molecules[find(PBS.parent, i)];
            if (molecule.atomCount == 1) {  
                const char* elem = atoms[i].element.c_str();
                drawCenteredLabel(atoms[i].pos.x, atoms[i].pos.y, elem, settings, RELATIVE_TEXT_HEIGHT_SOLO, Fade(GRAY, 0.8f));
            } else {
                std::string form = molecule.formula;
                std::string label = "[ " + form + "]";
                if (molecule.isProtein) label += " - Prtn.";
                const char* moleculeLabel = label.c_str();
                molecule.labelDrawn = true;
                drawCenteredLabel(molecule.centeroid.x, molecule.maxY, moleculeLabel, settings, RELATIVE_TEXT_HEIGHT_MOLECULE, Fade(GRAY, 0.8f));
            }
        }

        GuiToggleGroup({5, 565, 80, 30}, "Hydrogen;Carbon;Nitrogen;Oxygen;Phosphorus;Sulfur", &PBS.activeElement);
        bool main = GuiButton({10, 10, 30, 30}, "<");
        if (!PBS.isSaving) PBS.isSaving = GuiButton({715, 10, 75, 25}, "Save");
        else {
            GuiTextBox({700, 10, 90, 25}, PBS.saveBuf, settings.fontSize, true);
            bool save = GuiButton({700, 45, 90, 25}, "Save");
            if (save && PBS.saveBuf && PBS.saveBuf[0] != '\0') {
                std::filesystem::path home = getUserHomeDir();
                std::filesystem::path fpath = home / "Simulator";

                std::filesystem::create_directories(fpath);

                std::string fname = std::string(PBS.saveBuf) + ".protein";
                std::filesystem::path finalPath = fpath / fname;
                std::ofstream proteinStream(finalPath);

                if (!proteinStream.is_open()) {
                    std::cerr << "Failed to initialize protein write at protein file.\n";
                    exit(EXIT_FAILURE);
                }

                proteinStream << PBS.saveBuf << "\n";
                proteinStream << atoms.size() << "\n"; 
                for (size_t i = 0; i < atoms.size(); i++) {
                    proteinStream
                        << atoms[i].protons << ","
                        << atoms[i].neutrons << ","
                        << atoms[i].pos.x << ","
                        << atoms[i].pos.y << "\n";
                }
                proteinStream << PBS.bonds.size() << "\n";
                for (auto& bond : PBS.bonds) {
                    proteinStream 
                        << bond.first << ","
                        << bond.second << "\n";
                }

                PBS.isSaving = false;
            }
        }

        EndDrawing();

        if (main) return MAIN;

        frame++;
    }

    return QUIT;
}

SimScreen runBlankAtomSimulator() {
    SetWindowTitle(BLANK_SIMULATOR_WINDOW_TITLE);

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

    if (DEBUG_CONSOLE == 1 || DEBUG_CONSOLE == 3) { std::cout << "SEED: " << seed << "\n"; }

    // Other miscellaneous initializations
    int frame = 0;
    bool showSettings = false;
    std::vector<int> parent(atoms.size());
    std::vector<std::pair<size_t, size_t>> toFuse;
    std::vector<bool> fusing(atoms.size(), false);

    while (!WindowShouldClose()) {
        std::chrono::high_resolution_clock::time_point loopStartProfiler;
        if (DEBUG_PROFILER) loopStartProfiler = std::chrono::high_resolution_clock::now();

        std::chrono::high_resolution_clock::time_point startResetLoopProfiler;
        if (DEBUG_PROFILER) startResetLoopProfiler = std::chrono::high_resolution_clock::now();
        // Reset the parent array
        for (size_t i = 0; i < parent.size(); i++) {
            parent[i] = i;
        }
        
        for (size_t i = 0; i < fusing.size(); i++) {
            fusing[i] = false;
        }
        if (DEBUG_PROFILER) {
            auto endResetLoopProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsResetLoopProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endResetLoopProfiler - startResetLoopProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: RESET LOOPS: " << microsecondsResetLoopProfiler << "\n";
        }

        std::chrono::high_resolution_clock::time_point startPosProfiler;
        if (DEBUG_PROFILER) startPosProfiler = std::chrono::high_resolution_clock::now();
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
        if (DEBUG_PROFILER) {
            auto endPosProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsPosProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endPosProfiler - startPosProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: POS: " << microsecondsPosProfiler << "\n";
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
            if (frame % 60 == 0) std::cout << "PROFILER: N2: " << microsecondsN2Profiler << "\n";
        }

        std::chrono::high_resolution_clock::time_point startFuseAppProfiler;
        if (DEBUG_PROFILER) startFuseAppProfiler = std::chrono::high_resolution_clock::now();
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

        if (DEBUG_PROFILER) {
            auto endFuseAppProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsFuseAppProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endFuseAppProfiler - startFuseAppProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: FUSE APP.: " << microsecondsFuseAppProfiler << "\n";
        }

        std::chrono::high_resolution_clock::time_point startRDFProfiler;
        if (DEBUG_PROFILER) startRDFProfiler = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i].refreshDerivedFields(elementTable);
        }
        if (DEBUG_PROFILER) {
            auto endRDFProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsRDFProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endRDFProfiler - startRDFProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: RDF: " << microsecondsRDFProfiler << "\n";
        }

        std::chrono::high_resolution_clock::time_point startFrictionDecayProfiler;
        if (DEBUG_PROFILER) startFrictionDecayProfiler = std::chrono::high_resolution_clock::now();
        // Apply friction and decay to each atom
        for (size_t i = 0; i < atoms.size(); i++) {
            atoms[i].vel = atoms[i].vel * settings.atomicFriction;
            applyDecay(rng, atoms[i], elementTable, settings);
        }
        if (DEBUG_PROFILER) {
            auto endFrictionDecayProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsFrictionDecayProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endFrictionDecayProfiler - startFrictionDecayProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: FRIC/DECAY: " << microsecondsFrictionDecayProfiler << "\n";
        }
        
        std::chrono::high_resolution_clock::time_point startBuildMolecProfiler;
        if (DEBUG_PROFILER) startBuildMolecProfiler = std::chrono::high_resolution_clock::now();
        auto molecules = buildMolecules(atoms, parent);
        if (DEBUG_PROFILER) {
            auto endBuildMolecProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsBuildMolecProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endBuildMolecProfiler - startBuildMolecProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: BUILD MOLEC.: " << microsecondsBuildMolecProfiler << "\n";
        }
        
        // Log to a CSV and Console if needed
        for (size_t i = 0; i < atoms.size(); i++) {
            if (DEBUG_CONSOLE == 1 || DEBUG_CONSOLE == 3) { 
                std::cout 
                    << "Atom " << atoms[i].id << " - " << atoms[i].element
                    << ": X, Y, Z = " << atoms[i].pos.x << ", " << atoms[i].pos.y << ", " << atoms[i].pos.z << 
                    " - VX, VY, VZ = " << atoms[i].vel.x << ", " << atoms[i].vel.y << ", " << atoms[i].vel.z 
                    << " - Mass: " << atoms[i].mass
                    << "\n";
            }

            if (DEBUG_CONSOLE == 2 || DEBUG_CONSOLE == 3) {
                std::cout
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
                std::cout
                    << "Molecule (root " << entry.first << "): " << molecule.formula
                    << " - Mass: " << molecule.mass
                    << "\n";
            }
        }

        std::chrono::high_resolution_clock::time_point startDrawProfiler;
        if (DEBUG_PROFILER) startDrawProfiler = std::chrono::high_resolution_clock::now();
        // Draw to the screen (TODO)
        BeginDrawing();
        ClearBackground(BLACK);
        
        for (size_t i = 0; i < atoms.size(); i++) {
            renderAtom(atoms[i], settings, frame);
        }

        for (size_t i = 0; i < atoms.size(); i++) {
            // Add labels to each element
            Molecule& molecule = molecules[find(parent, i)];
            if (molecule.atomCount == 1) {
                const char* sym = atoms[i].element.c_str();
                drawCenteredLabel(atoms[i].pos.x, atoms[i].pos.y, sym, settings, RELATIVE_TEXT_HEIGHT_SOLO, Fade(GRAY, 0.8f));
            } else if (!molecule.labelDrawn) {
                std::string form = molecule.formula;
                std::string label = "[ " + form + "]";
                if (molecule.isProtein) label += " - Prtn.";
                const char* moleculeLabel = label.c_str();
                molecule.labelDrawn = true;
                drawCenteredLabel(molecule.centeroid.x, molecule.centeroid.y, moleculeLabel, settings, RELATIVE_TEXT_HEIGHT_MOLECULE, Fade(GRAY, 0.8f));
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
            if (settings.simpleAtoms) GuiDisable();
            GuiSlider(
                {150, 280, 300, 20}, 
                "Electron Orbit Speed", 
                TextFormat("%.3f", settings.electronOrbitSpeed), 
                &settings.electronOrbitSpeed, 0.0f, 0.1f
            );
            if (settings.simpleAtoms) GuiEnable();
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
            GuiCheckBox(
                {150, 370, 20, 20},
                "Simple Atoms (performance boost)",
                &settings.simpleAtoms
            );
        }

        if (GuiButton({710, 560, 70, 20}, "Settings")) {
            showSettings = !showSettings;
        }
        
        bool menu = GuiButton({10, 10, 30, 30}, "<");

        EndDrawing();

        if (menu) return MAIN;

        if (DEBUG_PROFILER) {
            auto endDrawProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsDrawProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endDrawProfiler - startDrawProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: DRAW: " << microsecondsDrawProfiler << "\n";
        }

        if (DEBUG_PROFILER) {
            auto endLoopProfiler = std::chrono::high_resolution_clock::now();
            auto microsecondsLoopProfiler = std::chrono::duration_cast<std::chrono::microseconds>(endLoopProfiler - loopStartProfiler).count();
            if (frame % 60 == 0) std::cout << "PROFILER: LOOP: " << microsecondsLoopProfiler << "\n";
        }

        frame++;
    }
    
    return QUIT;
}