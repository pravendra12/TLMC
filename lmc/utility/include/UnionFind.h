// UnionFind.hpp
#ifndef LMC_UTILITY_INCLUDE_UNIONFIND_H_
#define LMC_UTILITY_INCLUDE_UNIONFIND_H_

#include <unordered_map>
using namespace std;

template <typename T>
class UnionFind
{
public:
    T find(const T &x);
    void union_sets(const T &x, const T &y);
    unordered_map<T, T> getEquivalenceMap();

private:
    unordered_map<T, T> parent;
};

#include "UnionFind.inl"

#endif
