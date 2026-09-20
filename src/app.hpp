#ifndef SRC_APP_HPP
#define SRC_APP_HPP

struct Settings {
    float atomicFriction = 0.975f;
    float forceStrength = 0.2f;
    int forceMinDist = 15;
    int forceMaxDist = 50;
    int bondDist = 22;
    float decayChance = 0.04f;
    int decayProtonThreshold = 84;
    int fusionProtonMax = 70;
    float temp = 1.0f;
    float electronOrbitSpeed = 0.015f;
    int fontSize = 8;
};

#endif