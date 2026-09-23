#ifndef SIM_RENDER_HPP
#define SIM_RENDER_HPP

#include "raylib.h"

#include "atom.hpp"
#include "../core/settings.hpp"

void renderAtom(const Atom& atom, const Settings& settings, int frame);
void drawCenteredLabel(float x, float y, const char* text, const Settings& settings, int yOffset, Color color);

#endif