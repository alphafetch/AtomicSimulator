#ifndef CORE_SETTINGS_HPP
#define CORE_SETTINGS_HPP

struct Settings {
    float atomicFriction = 0.975f;
    float forceStrength = 0.2f;
    float forceMinDist = 15;
    float forceMaxDist = 50;
    float bondDist = 22;
    float decayChance = 0.04f;
    int decayProtonThreshold = 84;
    int fusionProtonMax = 70;
    float temp = 1.0f;
    float electronOrbitSpeed = 0.015f;
    int fontSize = 8;
};

#endif
