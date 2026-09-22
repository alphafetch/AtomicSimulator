#include "render.hpp"

#include <numbers>

#include "raylib.h"

#include "../core/constants.hpp"
#include "../core/vector.hpp"
#include "utils.hpp"

void renderAtom(const Atom& atom, const Settings& settings, int frame) {
    if (settings.simpleAtoms) {
        DrawCircleV(
            {atom.pos.x, atom.pos.y},
            ELECTRON_FLOAT_RADIUS + ELECTRON_RENDER_RADIUS,
            GREEN
        );
    } else {
        DrawCircleLinesV(Vector2({atom.pos.x, atom.pos.y}), ELECTRON_FLOAT_RADIUS, Fade(RED, 0.5f));

        for (int j = 0; j < atom.subatomTier[2]; j++) {
            float theta = j * ((2 * std::numbers::pi) / atom.subatomTier[2]) + frame * settings.electronOrbitSpeed;
            
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
}