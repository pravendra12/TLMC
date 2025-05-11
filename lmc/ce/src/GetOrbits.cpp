#include "GetOrbits.h"
#include "PrintUtility.h"

map<string, vector<vector<size_t>>> GetOrbits(
    const Config &config,
    const size_t &maxClusterSize,
    const size_t &maxBondOrder,
    const vector<vector<size_t>> &equivalentLatticeIdsVector)
{
  // Map lattice Id to equivalent site index
  /*
  unordered_map<size_t, size_t> latticeIdToIndexMap;
  for (size_t i = 0; i < symmetricallSortedVector.size(); i++)
  {
    latticeIdToIndexMap.insert(make_pair(symmetricallSortedVector[i], i));
  }
  */

  unordered_map<size_t, size_t> latticeIdToOrbitMap;

  for (size_t idx = 0; idx < equivalentLatticeIdsVector.size(); idx++)
  {
    for (auto latticeId : equivalentLatticeIdsVector[idx])
    {
      latticeIdToOrbitMap.insert(make_pair(latticeId, idx));
    }
  }

  vector<size_t> allLatticeIdVector;
  for (const auto &group : equivalentLatticeIdsVector)
  {
    allLatticeIdVector.insert(allLatticeIdVector.end(), group.begin(), group.end());
  }

  auto latticeClusterHashSet = FindClustersWithinAllowedSites(
      config,
      maxClusterSize,
      maxBondOrder,
      allLatticeIdVector);

  map<string, vector<vector<size_t>>> orbitMap;
  // Iterate over all the clusters

  for (const auto &latticeCluster : latticeClusterHashSet)
  {

    auto latticeIdVector = latticeCluster.GetLatticeIdVector();

    // vector<size_t> indexVector;
    vector<size_t> orbitVector;

    // corresopind index vector and orbit

    int orbitSum = 0;
    for (auto id : latticeIdVector)
    {
      // auto index = latticeIdToIndexMap.at(id);
      auto orbit = latticeIdToOrbitMap.at(id);

      // indexVector.emplace_back(index);
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