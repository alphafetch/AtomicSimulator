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
#include "../core/settings.hpp"

#include "../sim/screens.hpp"
#include "../sim/atom.hpp"
#include "../sim/utils.hpp"
#include "../sim/bondFind.hpp"
#include "../sim/molecule.hpp"
#include "../sim/render.hpp"
#include "../sim/proteinBuilder.hpp"

using std::cout;

int main() {
    InitWindow(WINDOW_WIDTH, WINDOW_HEIGHT, MAIN_WINDOW_TITLE);
    if (!DEBUG_PROFILER) { SetTargetFPS(60); }
    else { SetTargetFPS(12000); }

    SimScreen currentScreen = MAIN;
    do {
        if (currentScreen == MAIN) {
            currentScreen = runMainMenu();
        } else if (currentScreen == ATOM_SIMULATION) {
            currentScreen = runAtomSimulator();
        } else if (currentScreen == PROTEIN_BUILDER) {
            currentScreen = runProteinBuilder();
        }
    } while (currentScreen != QUIT);

    CloseWindow();
    exit(EXIT_SUCCESS);
}
