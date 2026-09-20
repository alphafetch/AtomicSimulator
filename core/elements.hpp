#ifndef CORE_ELEMENTS_HPP
#define CORE_ELEMENTS_HPP

#include <string>
#include <unordered_map>

struct Element {
    std::string symbol;
    std::string name;
    int neutrons;      // standard (most common isotope) neutron count
    long abundancePpm; // real cosmic (mass fraction) abundance, in ppm;
                       // elements without precise cited data use a small
                       // floor value (1) rather than 0, so they stay possible
};

// Loads a CSV with columns: atomic_number,symbol,name,protons,neutrons,abundance_ppm
// and returns a lookup from proton count -> Element. Proton count is what
// actually determines an element's identity (that's the definition of
// atomic number) - neutron/electron counts are ignored for matching since
// they vary by isotope/ion without changing what element something is.
std::unordered_map<int, Element> loadElementTable(const std::string& path);

#endif
