/*
 * sorters.cpp
 *
 *  Created on: Dec 31, 2012
 *      Author: snikolenko
 */
#include "new_electric.hpp"
#include "sorters.hpp"

bool Sorter::pairSecond(const pair<size_t, float> & p1, const pair<size_t, float> & p2) { return p1.second < p2.second; }

bool Sorter::demandsLength(const User & u1, const User & u2) { return u1.l > u2.l; }
bool Sorter::demandsAvailLength(const User & u1, const User & u2) { return u1.s_avail < u2.s_avail; }
bool Sorter::demandsRatio(const User & u1, const User & u2) { return u1.l / (float)u1.s_avail > u2.l / (float)u2.s_avail; }

bool Sorter::avgLoads(const pair<size_t, float> & u1, const pair<size_t, float> & u2) { return u1.second > u2.second; }
bool Sorter::avgLoadsReverse(const pair<size_t, float> & u1, const pair<size_t, float> & u2) { return u1.second < u2.second; }
