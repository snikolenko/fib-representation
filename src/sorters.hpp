/*
 * sorters.hpp
 *
 *  Created on: Dec 31, 2012
 *      Author: snikolenko
 */

#ifndef SORTERS_HPP_
#define SORTERS_HPP_

using namespace std;

class User;

class Sorter {
public:
	static bool pairSecond(const pair<size_t, float> & p1, const pair<size_t, float> & p2);

	static bool demandsLength(const User & u1, const User & u2);
	static bool demandsAvailLength(const User & u1, const User & u2);
	static bool demandsRatio(const User & u1, const User & u2);

	static bool avgLoads(const pair<size_t, float> & u1, const pair<size_t, float> & u2);
	static bool avgLoadsReverse(const pair<size_t, float> & u1, const pair<size_t, float> & u2);
};

class DemandSorter {
public:
	virtual void sort( vector<User> & v ) = 0;
};


class MLDFSorter : public DemandSorter {
public:
	virtual void sort( vector<User> & v ) { std::sort(v.begin(), v.end(), Sorter::demandsLength); }
};

class MRDFSorter : public DemandSorter {
public:
	virtual void sort( vector<User> & v ) { std::sort(v.begin(), v.end(), Sorter::demandsAvailLength); }
};

class MRatioDFSorter : public DemandSorter {
public:
	virtual void sort( vector<User> & v ) { std::sort(v.begin(), v.end(), Sorter::demandsRatio); }
};

#endif /* SORTERS_HPP_ */
