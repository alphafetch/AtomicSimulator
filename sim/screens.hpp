#ifndef SIM_SCREENS_HPP
#define SIM_SCREENS_HPP

typedef enum SimScreen { 
    MAIN, BLANK_SIMULATION, ATOM_SIMULATION, PROTEIN_BUILDER, QUIT
} SimScreen;

SimScreen runMainMenu();
SimScreen runAtomSimulator();
SimScreen runProteinBuilder();
SimScreen runBlankAtomSimulator();

#endif