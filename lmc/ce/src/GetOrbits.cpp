#include "GetOrbits.h"
#include "PrintUtility.h"

map<string, vector<vector<size_t>>> GetOrbits(
    const Config &config,
    const size_t &maxClusterSize,
    const size_t &maxBondOrder,
    const vector<vector<size_t>> &equivalentSiteEncoding,
    const vector<size_t> &symmetricallSortedVector)
{
  // Map lattice Id to equivalent site index

  unordered_map<size_t, size_t> latticeIdToIndexMap;
  for (size_t i = 0; i < symmetricallSortedVector.size(); i++)
  {
    latticeIdToIndexMap.insert(make_pair(symmetricallSortedVector[i], i));
  }

  unordered_map<size_t, size_t> indexToOrbitMap;

  for (size_t i = 0; i < equivalentSiteEncoding.size(); i++)
  {
    for (auto j : equivalentSiteEncoding[i])
    {
      indexToOrbitMap.insert(make_pair(j, i));
    }
  }

  auto latticeClusterHashSet = FindClustersWithinAllowedSites(
      config,
      maxClusterSize,
      maxBondOrder,
      symmetricallSortedVector);

  map<string, vector<vector<size_t>>> orbitMap;
  // Iterate over all the clusters

  for (const auto &latticeCluster : latticeClusterHashSet)
  {
    
    auto latticeIdVector = latticeCluster.GetLatticeIdVector();

    vector<size_t> indexVector;
    vector<size_t> orbitVector;

    // corresopind index vector and orbit

    int orbitSum = 0;
    for (auto id : latticeIdVector)
    {
      auto index = latticeIdToIndexMap.at(id);
      auto orbit = indexToOrbitMap.at(index);

      indexVector.emplace_back(index);
      orbitVector.emplace_back(orbit);

      orbitSum += orbit;
    }
    // print1DVector(indexVector);
    // print1DVector(orbitVector);

    // Sort orbits to create a unique key
    sort(orbitVector.begin(), orbitVector.end());
    string orbitKey;
    for (auto orbit : orbitVector)
    {
      orbitKey += to_string(orbit);
    }

    // cout << orbitKey << endl;

    // Handle empty clusters separately
    if (indexVector.empty())
    {
      orbitMap["-1"].emplace_back(indexVector);
    }
    else
    {
      orbitMap[orbitKey].emplace_back(indexVector);
    }
  }

  return orbitMap;
}