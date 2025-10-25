#include "B2OrderParameter.h"

B2OrderParameter::B2OrderParameter(const Config &config) : config(config) {
  InitializeAlphaLatticeSites();
  InitializeBetaLatticeSites();
}

double
B2OrderParameter::GetB2OrderParameter(const TiledSupercell &tiledSupercell,
                                      const Element &element) {
  // fractional occupancy of element at alpha and beta sites
  auto alphaOccupancy = GetAlphaSiteOccupancy(tiledSupercell, element);

  auto betaOccupancy = GetBetaSiteOccupancy(tiledSupercell, element);

  double b2OrderParameter =
      (alphaOccupancy - betaOccupancy) / (alphaOccupancy + betaOccupancy);

  return b2OrderParameter;
}

pair<double, double> B2OrderParameter::GetAlphaBetaSiteOccupancy(
    const TiledSupercell &tiledSupercell, const Element &element) {
  // Number of occupied alpha sites by element
  size_t numOccupiedAlphaSites = 0;
  size_t numOccupiedBetaSites = 0;

  size_t numAlphaSites = alphaLatticeSites.size();
  size_t numBetaSites = betaLatticeSites.size();

  for (size_t smallConfigIdx = 0;
       smallConfigIdx < tiledSupercell.GetNumOfSmallConfig();
       ++smallConfigIdx) {
    for (const auto &alphaLatticeId : alphaLatticeSites) {
      auto siteElement = tiledSupercell.GetElementAtSite(
          LatticeSiteMapping(alphaLatticeId, smallConfigIdx));

      if (siteElement == element) {
        numOccupiedAlphaSites += 1;
      }
    }

    for (const auto &betaLatticeId : betaLatticeSites) {
      auto siteElement = tiledSupercell.GetElementAtSite(
          LatticeSiteMapping(betaLatticeId, smallConfigIdx));

      if (siteElement == element) {
        numOccupiedBetaSites += 1;
      }
    }
  }

  double alphaSiteOccupancy =
      double(numOccupiedAlphaSites) / double(numAlphaSites);
  double betaSiteOccupancy =
      double(numOccupiedBetaSites) / double(numBetaSites);

  return {alphaSiteOccupancy, betaSiteOccupancy};
}

double B2OrderParameter::GetBetaSiteOccupancy(const Element &element) {
  // Number of occupied beta sites by element
  size_t numOccupiedBSites = 0;
  size_t numBetaSites = betaLatticeSites.size();

  for (const auto &id : betaLatticeSites) {
    auto siteElement = config.GetElementOfLattice(id);

    if (siteElement == element) {
      numOccupiedBSites += 1;
    }
  }

  double betaSiteOccupancy = double(numOccupiedBSites) / double(numBetaSites);

  return betaSiteOccupancy;
}

unordered_set<size_t> B2OrderParameter::GetAlphaLatticeSites() {
  return alphaLatticeSites;
}

unordered_set<size_t> B2OrderParameter::GetBetaLatticeSites() {
  return betaLatticeSites;
}

bool B2OrderParameter::isB2Ordered(const Config &config, const size_t atomId) {
  // Element vacancy("X");

  // if (config.GetElementOfAtom(atomId) == vacancy)
  // {
  //   return false;
  // }

  Element centerElement = config.GetElementOfAtom(atomId);
  auto firstNN = config.GetNeighborAtomIdVectorOfAtom(atomId, 1);

  Element neighborElement = config.GetElementOfAtom(firstNN[0]);

  for (size_t i = 0; i < firstNN.size(); ++i) {
    if (config.GetElementOfAtom(firstNN[i]) != neighborElement) {
      return false;
    }
  }

  if (neighborElement == centerElement)
    return false;

  // second NN should all be equal to centerElement for perfect B2

  auto secondNN = config.GetNeighborAtomIdVectorOfAtom(atomId, 2);

  // case of defected B2

  for (size_t i = 0; i < secondNN.size(); ++i) {
    if (config.GetElementOfAtom(secondNN[i]) != centerElement) {
      return false;
    }
  }

  return true;
}

void B2OrderParameter::InitializeAlphaLatticeSites() {
  // Precompute first nearest neighbors map for constant-time lookups
  const auto &firstNNList = config.GetNeighborLists()[0];
  std::vector<std::unordered_set<size_t>> firstNNMap(firstNNList.size());

  for (size_t i = 0; i < firstNNList.size(); ++i) {
    firstNNMap[i] = std::unordered_set<size_t>(firstNNList[i].begin(),
                                               firstNNList[i].end());
  }

  auto secondNNList = config.GetNeighborLists()[1];

  // To avoid processing same pair twice
  std::set<std::pair<size_t, size_t>> visited_pairs;

  for (size_t id1 = 0; id1 < secondNNList.size(); ++id1) {
    const auto &secondNN = secondNNList[id1];

    for (size_t id2 : secondNN) {
      // Skip if this pair was already processed
      auto pair_key = std::minmax(id1, id2);
      if (visited_pairs.count(pair_key) > 0)
        continue;
      visited_pairs.insert(pair_key);

      bool id1Valid = true;
      bool id2Valid = true;

      // Check if id1 has any first NN in alphaLatticeSites
      for (size_t id3 : alphaLatticeSites) {
        if (firstNNMap[id1].count(id3)) {
          id1Valid = false;
          break;
        }
      }

      // Check if id2 has any first NN in alphaLatticeSites
      for (size_t id3 : alphaLatticeSites) {
        if (firstNNMap[id2].count(id3)) {
          id2Valid = false;
          break;
        }
      }

      if (id1Valid) {
        alphaLatticeSites.emplace(id1);
      }
      if (id2Valid) {
        alphaLatticeSites.emplace(id2);
      }
    }
  }
}

void B2OrderParameter::InitializeBetaLatticeSites() {
  // betaSites = allSites - alphaSites
  for (size_t id = 0; id < numLattice; id++) {
    // O(1) Time Complexity
    if (alphaLatticeSites.find(id) == alphaLatticeSites.end()) {
      betaLatticeSites.emplace(id);
    }
  }
}
