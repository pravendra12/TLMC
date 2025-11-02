#include "SubLatticeOccupancy.h"

SubLatticeOccupancy::SubLatticeOccupancy(
    const TiledSupercell &tiledSupercell) : alphaAtomIds_(InitializeAlphaAtomIds(tiledSupercell))
{
}

map<Element, double> SubLatticeOccupancy::ComputeOrderParameter(
    const map<Element, pair<double, double>> &elementOccupancies) const
{
  map<Element, double> elementOrderParameter;

  for (const auto &[element, occupancies] : elementOccupancies)
  {
    const double alphaOccupancy = occupancies.first;
    const double betaOccupancy = occupancies.second;

    double eta = 0.0;
    const double denom = alphaOccupancy + betaOccupancy;

    if (denom > 1e-12) // avoid division by zero
      eta = (alphaOccupancy - betaOccupancy) / denom;

    elementOrderParameter[element] = eta;
  }

  return elementOrderParameter;
}


map<Element, pair<double, double>> SubLatticeOccupancy::GetAlphaBetaSiteOccupancy(
    const vector<size_t> &atomIndicesVector) const
{
  const size_t totalNumSites = atomIndicesVector.size();

  const size_t numAlphaSites = totalNumSites / 2;
  const size_t numBetaSites = numAlphaSites;

  // Counters for all elements: element_id -> (num_alpha, num_beta)
  map<Element, pair<size_t, size_t>> elementCounts;

  for (size_t atomId = 0; atomId < totalNumSites; ++atomId)
  {
    const size_t atomicNumber = atomIndicesVector[atomId];
    const bool isAlpha = (alphaAtomIds_.find(atomId) != alphaAtomIds_.end());

    auto &counts = elementCounts[Element(atomicNumber)];
    if (isAlpha)
      ++counts.first; // alpha count
    else
      ++counts.second; // beta count
  }

  // Convert counts → occupancy fractions
  map<Element, pair<double, double>> elementOccupancies;

  for (const auto &[element, counts] : elementCounts)
  {
    double alphaOcc = static_cast<double>(counts.first) / static_cast<double>(numAlphaSites);
    double betaOcc = static_cast<double>(counts.second) / static_cast<double>(numBetaSites);
    elementOccupancies[element] = {alphaOcc, betaOcc};
  }

  return elementOccupancies;
}

const unordered_set<size_t> &SubLatticeOccupancy::GetAlphaAtomIds() const
{
  return alphaAtomIds_;
}

unordered_set<size_t> SubLatticeOccupancy::InitializeAlphaAtomIds(
    const TiledSupercell &tiledSupercell)
{
  // Use small config to get the alpha sites then assign the cubeId later
  const Config &smallConfig = tiledSupercell.GetSmallConfig();

  unordered_set<size_t> alphaLatticeSites;

  // Precompute first nearest neighbors map for constant-time lookups
  const auto &firstNNList = smallConfig.GetNeighborLists()[0];
  vector<unordered_set<size_t>> firstNNMap(firstNNList.size());

  for (size_t i = 0; i < firstNNList.size(); ++i)
  {
    firstNNMap[i] = unordered_set<size_t>(firstNNList[i].begin(),
                                          firstNNList[i].end());
  }

  auto secondNNList = smallConfig.GetNeighborLists()[1];

  // To avoid processing same pair twice
  set<pair<size_t, size_t>> visited_pairs;

  for (size_t id1 = 0; id1 < secondNNList.size(); ++id1)
  {
    const auto &secondNN = secondNNList[id1];

    for (size_t id2 : secondNN)
    {
      // Skip if this pair was already processed
      auto pair_key = minmax(id1, id2);
      if (visited_pairs.count(pair_key) > 0)
        continue;
      visited_pairs.insert(pair_key);

      bool id1Valid = true;
      bool id2Valid = true;

      // Check if id1 has any first NN in alphaLatticeSites
      for (size_t id3 : alphaLatticeSites)
      {
        if (firstNNMap[id1].count(id3))
        {
          id1Valid = false;
          break;
        }
      }

      // Check if id2 has any first NN in alphaLatticeSites
      for (size_t id3 : alphaLatticeSites)
      {
        if (firstNNMap[id2].count(id3))
        {
          id2Valid = false;
          break;
        }
      }

      if (id1Valid)
      {
        alphaLatticeSites.emplace(id1);
      }
      if (id2Valid)
      {
        alphaLatticeSites.emplace(id2);
      }
    }
  }

  const size_t numOfSmallConfigs = tiledSupercell.GetNumOfSmallConfig();

  unordered_set<size_t> alphaAtomIds;

  for (size_t cId = 0; cId < numOfSmallConfigs; cId++)
  {
    for (const size_t latticeId : alphaLatticeSites)
    {

      auto atomId = tiledSupercell.GetAtomIdFromLatticeAndConfigId(
          LatticeSiteMapping(
              latticeId,
              cId));

      alphaAtomIds.insert(atomId);
    }
  }

  return alphaAtomIds;
}
