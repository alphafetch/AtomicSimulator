#ifndef CORE_ELEMENTS_HPP
#define CORE_ELEMENTS_HPP

#include <string>
#include <unordered_map>

struct Element {
    std::string symbol;
    std::string name;
    int neutrons; 
    long abundancePpm; 
};

std::unordered_map<int, Element> loadElementTable(const std::string& path);

#endif
