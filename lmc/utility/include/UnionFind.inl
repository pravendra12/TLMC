#include "UnionFind.h"
template <typename T>
T UnionFind<T>::find(const T& x) {
    if (parent.find(x) == parent.end())
        parent[x] = x;
    if (parent[x] != x)
        parent[x] = find(parent[x]);
    return parent[x];
}

template <typename T>
void UnionFind<T>::union_sets(const T& x, const T& y) {
    T rootX = find(x);
    T rootY = find(y);
    parent[rootX] = rootY;
}

template <typename T>
inline unordered_map<T, T> UnionFind<T>::getEquivalenceMap()
{
  return parent;
}
