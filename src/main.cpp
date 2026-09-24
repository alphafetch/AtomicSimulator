#include <cstdlib>

#include "raylib.h"
#define RAYGUI_IMPLEMENTATION
#include "../external/raygui.h"

#include "../core/constants.hpp"

#include "../sim/screens.hpp"

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
