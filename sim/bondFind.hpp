#ifndef SIM_UNIONFIND_HPP
#define SIM_UNIONFIND_HPP

#include <vector>

int find(std::vector<int>& parent, int x);

void bond(std::vector<int>& parent, int xid, int yid);

#endif