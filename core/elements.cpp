#include "elements.hpp"

#include <fstream>
#include <sstream>
#include <iostream>

std::unordered_map<int, Element> loadElementTable(const std::string& path) {
    std::unordered_map<int, Element> table;

    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Failed to open element table at " << path << "\n";
        return table;
    }

    std::string line;
    std::getline(file, line); // skip header row ("atomic_number,symbol,name,protons,neutrons")

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string field;

        std::getline(ss, field, ',');   // atomic_number column
        int protons = std::stoi(field);

        Element el;
        std::getline(ss, el.symbol, ',');
        std::getline(ss, el.name, ',');

        std::getline(ss, field, ','); // protons column - redundant with atomic_number, unused
        std::getline(ss, field, ','); // neutrons column
        el.neutrons = std::stoi(field);

        table[protons] = el;
    }

    return table;
}
