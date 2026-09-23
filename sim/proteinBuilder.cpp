#include "proteinBuilder.hpp"

void ProteinBuilderStorage::resizeParent(size_t s) { 
    size_t old = this->parent.size();
    this->parent.resize(s);
    for (size_t i = old; i < s; i++) {
        this->parent[i] = i;
    }
}