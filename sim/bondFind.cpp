#include "bondFind.hpp"

int find(std::vector<int>& parent, int x) {
    while (parent[x] != x) {
        x = parent[x];
    }

    return x;
}

void bond(std::vector<int>& parent, int xid, int yid) {
    int rootX = find(parent, xid);
    int rootY = find(parent, yid);

    if (rootX != rootY) {
        parent[rootY] = rootX;
    }
}