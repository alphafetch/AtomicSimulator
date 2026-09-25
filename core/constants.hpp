#ifndef CORE_CONSTANTS_HPP
#define CORE_CONSTANTS_HPP

#define COUNT 300

#define NEUTRON_MASS_KG 1.67492749804e-27
#define ELECTRON_MASS_KG 9.1093837e-31
#define PROTON_MASS_KG 1.67262192e-27
#define SUBATOMIC_RENDER_RADIUS 3.0f
#define NUCLEAR_FLOAT_RADIUS 3
#define ELECTRON_RENDER_RADIUS 2
#define ELECTRON_FLOAT_RADIUS 9
#define SELECT_RING_ADDITIONAL_RAD 2
#define BOND_THICKNESS 2
#define MIN_PROTEIN_SIZE 6

#define WINDOW_WIDTH 800
#define WINDOW_HEIGHT 600
#define MAIN_WINDOW_TITLE "[ Simulator ]"
#define ATOM_SIMULATOR_WINDOW_TITLE "[ Atomic Simulator ]"
#define PROTEIN_BUILDER_WINDOW_TITLE "[ Protein Builder ]"
#define BLANK_SIMULATOR_WINDOW_TITLE "[ Blank Simulation ]"

#define RELATIVE_TEXT_HEIGHT_SOLO -20
#define RELATIVE_TEXT_HEIGHT_MOLECULE -35

// DEBUG_CONSOLE: 0 - None
//                1 - Log atomic positions
//                2 - Log bonds
//                3 - All logs
//                4 - Molecule logging
#define DEBUG_CONSOLE 0
#define DEBUG_CSV 0
#define DEBUG_PROFILER 0

#endif
