#include "GetOrbits.h"
#include "PrintUtility.h"

#include <chrono>

using namespace chrono;

map<string, vector<vector<size_t>>> GetOrbits(
    const Config &config,
    const size_t &maxClusterSize,
    const size_t &maxBondOrderOfClusters,
    const vector<size_t> &ssVector,
    const vector<vector<size_t>> &equivalentSitesEncoding)
{
  // Map lattice Id to equivalent site index
  unordered_map<size_t, size_t> latticeIdToIndexMap;
  latticeIdToIndexMap.reserve(ssVector.size());

  for (size_t i = 0; i < ssVector.size(); i++)
  {
    latticeIdToIndexMap.insert(make_pair(ssVector[i], i));
  }

  unordered_map<size_t, size_t> indexToOrbitMap;
  indexToOrbitMap.reserve(ssVector.size());

  for (size_t idx = 0; idx < equivalentSitesEncoding.size(); idx++)
  {
    for (auto eqSiteIndex : equivalentSitesEncoding[idx])
    {
      indexToOrbitMap.insert(make_pair(eqSiteIndex, idx));
    }
  }

  auto start = high_resolution_clock::now();

  auto latticeClusterHashSet = FindClustersWithinAllowedSites(
      config,
      maxClusterSize,
      maxBondOrderOfClusters,
      ssVector);

  auto end = high_resolution_clock::now();

  auto duration = duration_cast<microseconds>(end - start);

  // cout << "Time to compute FindClustersWithinAllowedSites: " << duration.count() << " microseconds" << endl;

  // Use of map for consistent order
  map<string, vector<vector<size_t>>> orbitMap;

  // Iterate over all the lattice clusters
  for (const auto &latticeCluster : latticeClusterHashSet)
  {
    const auto &latticeIdVector = latticeCluster.GetLatticeIdVector();

    vector<size_t> orbitVector;
    orbitVector.reserve(latticeIdVector.size());

    for (const auto &latticeId : latticeIdVector)
    {
      const auto &index = latticeIdToIndexMap.at(latticeId);
      const auto &orbit = indexToOrbitMap.at(index);
      orbitVector.emplace_back(orbit);
    }

    sort(orbitVector.begin(), orbitVector.end());

    // Unique key for the orbit based on the order of the equivalent sites encoding
    string orbitKey;
    orbitKey.reserve(orbitVector.size()); 
    for (const auto &orbit : orbitVector)
      orbitKey += to_string(orbit);

    // Store the orbits
    if (latticeIdVector.empty())
    {
      orbitMap["-1"].emplace_back(latticeIdVector);
    }
    else
    {
      orbitMap[orbitKey].emplace_back(latticeIdVector);
    }
  }

  return orbitMap;
}