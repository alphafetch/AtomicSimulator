#include "screens.hpp"

#include <cstdlib>

#include "raylib.h"
#include "../external/raygui.h"

SimScreen MAINMenu() {
    while (!WindowShouldClose()) {
        BeginDrawing();

        ClearBackground(DARKGRAY);

        if (GuiButton({340, 265, 120, 40}, "Atom Simulation")) {
            return ATOM_SIMULATION;
        }

        if (GuiButton({340, 315, 120, 40}, "Protein Builder")) {
            return PROTEIN_BUILDER;
        }

        EndDrawing();
    }

    CloseWindow();
    exit(EXIT_SUCCESS);
}