#include "ShortRangeOrderTLMC.h"

ShortRangeOrderTLMC::ShortRangeOrderTLMC(
    const TiledSupercell &tiledSupercell)
    : tiledSupercell_(tiledSupercell),
      compositionMap_(tiledSupercell_.GetConcentrationMap())
{
}

// Compute Warren–Cowley SRO
map<string, double> ShortRangeOrderTLMC::ComputeWarrenCowley(const size_t shellNumber)
{
  // Get directed pair counts: center -> (neighbor -> count)
  PairCountByCenter pairCounts = GetElementPairCountMap(shellNumber);

  map<string, double> alphaMap;

  for (const auto &center_kv : pairCounts)
  {
    const string &elem_i = center_kv.first;   // center element
    const auto &neigh_map = center_kv.second; // neighbor -> count

    // total neighbors counted for centers of type i (sum over all neighbor types)
    size_t totalNeighbors_i = 0;
    for (const auto &kv : neigh_map)
      totalNeighbors_i += kv.second;

    if (totalNeighbors_i == 0)
    {
      // nothing to compute for this center element
      continue;
    }

    // compute probabilities and alpha for each neighbor element j
    for (const auto &kv : neigh_map)
    {
      const string &elem_j = kv.first;
      size_t count_ij = kv.second;

      double P_ij = static_cast<double>(count_ij) / static_cast<double>(totalNeighbors_i);

      // access composition map. Your compositionMap_ used Element as key earlier,
      // so we construct Element(elem_j). If your map uses strings, adapt here.
      double c_j = 0.0;
      try
      {
        c_j = compositionMap_.at(Element(elem_j));
      }
      catch (...)
      {
        // if Element(...) key not found, try fallback (if compositionMap_ uses strings)
        // Uncomment the following line if you maintain a string->double map instead:
        // c_j = compositionMapStr_.at(elem_j);
        throw std::runtime_error("Error in `ShortRangeOrderTLMC::ComputeWarrenCowley`: Composition for element " + elem_j + " not found in compositionMap_");
      }

      // Avoid division by zero (shouldn't happen if composition is valid)
      if (c_j <= 0.0)
      {
        throw std::runtime_error("Error in `ShortRangeOrderTLMC::ComputeWarrenCowley`: Invalid composition (zero) for element " + elem_j);
      }

      double alpha_ij = (P_ij - compositionMap_.at(Element(elem_j))) /
                        (static_cast<double>(elem_i == elem_j) - compositionMap_.at(Element(elem_j)));

      alphaMap[elem_i + "-" + elem_j] = alpha_ij;
    }
  }

  return alphaMap;
}

// Serial Version
/*
PairCountByCenter ShortRangeOrderTLMC::GetElementPairCountMap(const size_t shellNumber)
{
  const auto numOfSmallConfigs = tiledSupercell_.GetNumOfSmallConfig();

  PairCountByCenter totalPairCounts;

  for (size_t smallConfigId = 0; smallConfigId < numOfSmallConfigs; ++smallConfigId)
  {
    auto elementPairCountPerCenter = ComputeElementPairCountsForSmallConfig(
        smallConfigId,
        shellNumber);

    // Merge into total nested map
    for (const auto &center_kv : elementPairCountPerCenter)
    {
      const string &center = center_kv.first;
      const auto &neigh_map = center_kv.second;
      auto &target_neigh_map = totalPairCounts[center]; // creates if not exists

      for (const auto &nkv : neigh_map)
      {
        target_neigh_map[nkv.first] += nkv.second;
      }
    }
  }

  return totalPairCounts;
}
*/

// Return nested map: center_element -> (neighbor_element -> count)
PairCountByCenter ShortRangeOrderTLMC::GetElementPairCountMap(const size_t shellNumber)
{
  const auto numOfSmallConfigs = tiledSupercell_.GetNumOfSmallConfig();

  // Thread-local storage
  std::vector<PairCountByCenter> threadPairCounts;

#pragma omp parallel
  {
    PairCountByCenter localPairCounts;

#pragma omp for nowait
    for (size_t smallConfigId = 0; smallConfigId < numOfSmallConfigs; ++smallConfigId)
    {
      auto elementPairCountPerCenter = ComputeElementPairCountsForSmallConfig(
          smallConfigId,
          shellNumber);

      // Merge into thread-local map
      for (const auto &center_kv : elementPairCountPerCenter)
      {
        const string &center = center_kv.first;
        const auto &neigh_map = center_kv.second;
        auto &target_neigh_map = localPairCounts[center];

        for (const auto &nkv : neigh_map)
          target_neigh_map[nkv.first] += nkv.second;
      }
    }

#pragma omp critical
    {
      threadPairCounts.push_back(std::move(localPairCounts));
    }
  }

  // Merge all thread-local maps into total
  PairCountByCenter totalPairCounts;
  for (const auto &threadMap : threadPairCounts)
  {
    for (const auto &center_kv : threadMap)
    {
      const string &center = center_kv.first;
      const auto &neigh_map = center_kv.second;
      auto &target_neigh_map = totalPairCounts[center];

      for (const auto &nkv : neigh_map)
        target_neigh_map[nkv.first] += nkv.second;
    }
  }

  return totalPairCounts;
}

// Now ComputeElementPairCountsForSmallConfig returns directed counts:
// center_element -> (neighbor_element -> count)
PairCountByCenter ShortRangeOrderTLMC::ComputeElementPairCountsForSmallConfig(
    const size_t smallConfigId,
    const size_t shellNumber)
{
  const size_t numSitesInSmallConfig = tiledSupercell_.GetNumOfSitesPerSmallConfig();

  // neighbours mapping for small config
  const auto &neighboursOfSmallConfig = tiledSupercell_.GetCube().GetNeighbors(smallConfigId);

  PairCountByCenter elementPairMap;

  for (size_t latticeId = 0; latticeId < numSitesInSmallConfig; ++latticeId)
  {
    const LatticeSiteMapping latticeSiteId(latticeId, smallConfigId);
    auto centerElementObj = tiledSupercell_.GetElementAtSite(latticeSiteId);
    string centerElement = centerElementObj.GetElementString();

    if (centerElementObj == Element("X"))
      continue;

    auto nnEncodedVector = tiledSupercell_.GetNeighborLatticeIdVectorOfLattice(
        latticeId, shellNumber);

    for (const auto &nnEncodedId : nnEncodedVector)
    {
      size_t nnSmallConfigId = smallConfigId;
      if (nnEncodedId.encodedSmallConfigId != -1)
      {
        nnSmallConfigId = neighboursOfSmallConfig[size_t(nnEncodedId.encodedSmallConfigId)];
      }

      const LatticeSiteMapping nnLatticeSiteId(nnEncodedId.latticeId, nnSmallConfigId);
      auto neighbourElementObj = tiledSupercell_.GetElementAtSite(nnLatticeSiteId);

      if (neighbourElementObj == Element("X"))
        continue;

      string neighbourElement = neighbourElementObj.GetElementString();

      // increment directed count: center -> neighbour
      elementPairMap[centerElement][neighbourElement] += 1;
    }
  }

  return elementPairMap;
}

string ShortRangeOrderTLMC::GetElementPairString(
    const Element &element1,
    const Element &element2)
{
  if (element1 < element2)
  {
    return element1.GetElementString() + "-" + element2.GetElementString();
  }
  else
  {
    return element2.GetElementString() + "-" + element1.GetElementString();
  }
}
