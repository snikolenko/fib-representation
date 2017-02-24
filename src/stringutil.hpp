#pragma once

#include <algorithm>
#include <string>

#include <boost/algorithm/string.hpp>

using namespace std;

/************ FOR STRINGS *****************/

vector<string> & split(const string & s, char delim, vector<string> & elems) {
    stringstream ss(s);
    string item;
    elems.clear();
    if (s.size() <= 0) return elems;
    while(getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string & s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

// trim from start
static inline string &ltrim(string &s) {
        s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
        return s;
}

// trim from end
static inline string &rtrim(string &s) {
        s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
        return s;
}

// trim from both ends
static inline string &trim(string &s) {
        return ltrim(rtrim(s));
}

string join_vector(const std::vector<string> v, size_t start, size_t end) {
    string res = "";
    for (size_t i = start; i < end; ++i) {
        res += v[i] + ";";
    }
    return res;
}
