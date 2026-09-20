#ifndef CORE_ELEMENTS_HPP
#define CORE_ELEMENTS_HPP

#include <string>
#include <unordered_map>

struct Element {
    std::string symbol;
    std::string name;
    int neutrons; // standard (most common isotope) neutron count
};

// Loads a CSV with columns: atomic_number,symbol,name,protons,neutrons
// and returns a lookup from proton count -> Element. Proton count is what
// actually determines an element's identity (that's the definition of
// atomic number) - neutron/electron counts are ignored for matching since
// they vary by isotope/ion without changing what element something is.
std::unordered_map<int, Element> loadElementTable(const std::string& path);

#endif
