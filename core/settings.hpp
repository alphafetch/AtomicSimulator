#ifndef CORE_SETTINGS_HPP
#define CORE_SETTINGS_HPP

struct Settings {
    float atomicFriction = 0.975f;
    float forceStrength = 0.2f;
    float forceMinDist = 15.0f;
    float forceMaxDist = 50.0f;
    float bondDist = 22.0f;
    float decayChance = 0.04f;
    float decayProtonThreshold = 84.0f;
    float fusionProtonMax = 70.0f;
    float temp = 1.0f;
    float electronOrbitSpeed = 0.015f;
    float fontSize = 8.0f;
    bool  simpleAtoms = false;
};

#endif
