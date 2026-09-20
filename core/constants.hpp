#ifndef CORE_CONSTANTS_HPP
#define CORE_CONSTANTS_HPP

#define COUNT 200

#define NEUTRON_MASS_KG 1.67492749804e-27
#define ELECTRON_MASS_KG 9.1093837e-31
#define PROTON_MASS_KG 1.67262192e-27
#define SUBATOMIC_RENDER_RADIUS 3.0f
#define NUCLEAR_FLOAT_RADIUS 3
#define ELECTRON_RENDER_RADIUS 2
#define ELECTRON_FLOAT_RADIUS 9
#define ELECTRON_ORBIT_SPEED 0.015f
#define ATOMIC_FRICTION 0.975f
#define FORCE_STRENGTH 0.2f

#define FORCE_MIN_DIST 15
#define FORCE_MAX_DIST 50
#define BOND_DIST 22

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600
#define WINDOW_TITLE "Quantum Simulator"

#define RELATIVE_TEXT_HEIGHT_SOLO -20
#define RELATIVE_TEXT_HEIGHT_MOLECULE -35
#define FONT_SIZE 8

// DEBUG_CONSOLE: 0 - None
//                1 - Log atomic positions
//                2 - Log bonds
//                3 - All logs
//                4 - Molecule logging
#define DEBUG_CONSOLE 0
#define DEBUG_CSV 0

#endif
