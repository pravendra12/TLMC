/*******************************************************************************
 * Copyright (c) 2025. All rights reserved.
 * @Author: Pravendra Patel
 * @Date:    2025-06-01
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 *******************************************************************************/

/*! @file KRAPredictorTLMC.cpp
    @brief File contains implementation of Kinetically Resolved Activation Barrier
 */

#include "KRAPredictorTLMC.h"

KRAPredictorTLMC::KRAPredictorTLMC(
    const ClusterExpansionParameters &ceParams,
    const TiledSupercell &tiledSupercell) : maxBondOrder_(ceParams.GetMaxBondOrder("kra")),
                                            maxClusterSize_(
                                                ceParams.GetMaxClusterSize("ce")),
                                            referenceJumpDirection_(
                                                ceParams.GetReferenceJumpDirection()),
                                            atomicBasis_(
                                                ceParams.GetElementSet("kra"),
                                                ceParams.GetBasisType()),
                                            canonicalReferenceMap_(
                                                GetCenteredNeighborsAlongJumpDirection(
                                                    tiledSupercell.GetSmallConfig(),
                                                    maxBondOrder_,
                                                    referenceJumpDirection_)),
                                            symmetryOperations_(
                                                GetSymmetryOperations(
                                                    tiledSupercell.GetSmallConfig())),
                                            encodedOrbitsForPair_(
                                                GetLocalEncodedOrbitsForPair(
                                                    tiledSupercell.GetSmallConfig(),
                                                    maxBondOrder_,
                                                    maxClusterSize_,
                                                    referenceJumpDirection_,
                                                    true)),
                                            KECIsMap_(
                                                ceParams.GetKECIsMap())

{
  // Precompute the canonical sorted lattice Ids
  // And PairToIndexMap
  GetSymmetricallySortedLatticeIdsVectorMap(tiledSupercell);

  // Print KRAPredictorInfo
  PrintKRAPredictorInfo();
}

double KRAPredictorTLMC::GetKRA(
    const TiledSupercell &tiledSupercell,
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair) const
{
  const auto &elementFirst = tiledSupercell.GetElementAtSite(latticeSiteJumpPair.first);
  const auto &elementSecond = tiledSupercell.GetElementAtSite(latticeSiteJumpPair.second);

  // At least one element must be a vacancy
  if (elementFirst != Element("X") && elementSecond != Element("X"))
  {
    throw runtime_error("Error in `KRAPredictorTLMC::GetKRA`: Jump pair must contain at least one vacancy.");
  }

  Element migratingElement = (elementFirst == Element("X"))
                                 ? elementSecond
                                 : elementFirst;

  // Canonical Pair
  auto site1 = latticeSiteJumpPair.first;
  auto site2 = latticeSiteJumpPair.second;

  auto canonicalSiteJumpPair = (site1.latticeId < site2.latticeId)
                                   ? make_pair(site1, site2)
                                   : make_pair(site2, site1);

  pair<size_t, size_t> canonicalLatticeIdJumpPair(
      canonicalSiteJumpPair.first.latticeId,
      canonicalSiteJumpPair.second.latticeId);

  size_t indexOfPair = pairToIndexHashMap_.at(canonicalLatticeIdJumpPair);
  const auto &canonicalSortedLatticeSites = symmetricLatticeSiteEncodedMappings_[indexOfPair];

  // Now this is tricky part because here in the case of pair one also need to think
  // about which neighbour corresponds to which latticeId
  // if a id is neighbour of both then the atom at that latticeSiteMapping will be same
  // but the latticeSiteMappingEncoded will be different for both sites

  // This has been solved using NeighbourOfPair, now the canonicalSortedLatticeSites
  // contains the information about which latticeSite they belong to simple

  VectorXd correlationVector = GetCorrelationVector(
      tiledSupercell,
      canonicalSiteJumpPair,
      atomicBasis_,
      canonicalSortedLatticeSites,
      encodedOrbitsForPair_);

  // EKRA = J.Φ_α
  VectorXd kecis = KECIsMap_.at(migratingElement);

  if (kecis.size() != correlationVector.size())
  {
    throw runtime_error(
        "Error in `KRAPredictorTLMC::GetKRA`: migratingElement `" + migratingElement.GetElementString() +
        "` has kecis size (" + to_string(kecis.size()) +
        ") which does not match correlationVector size (" + to_string(correlationVector.size()) +
        ").");
  }

  double eKraValue = kecis.dot(correlationVector);

  return eKraValue;
}

void KRAPredictorTLMC::GetSymmetricallySortedLatticeIdsVectorMap(
    const TiledSupercell &tiledSupercell)
{
  const Config &smallConfig = tiledSupercell.GetSmallConfig();
  // only have the unique neighbours encoding in case of clash
  // keep any one as it does not matter as the element will be same if computed in
  // one way or other, since its first nn so this can be used
  const bool removeLatticeIdDuplicates = true;
  //

  // Number of lattice in the small config
  const size_t numLattices = smallConfig.GetNumLattices();

  auto numNeighboursPerSite = smallConfig.GetNeighborLatticeIdVectorOfLattice(0, 1).size();
  const size_t numUniquePairs = (numLattices * numNeighboursPerSite) / 2;

  // Declared as private members
  symmetricLatticeSiteEncodedMappings_.reserve(numUniquePairs);
  pairToIndexHashMap_.reserve(numUniquePairs);

  size_t numSortedSites = canonicalReferenceMap_.size();

  vector<NeighbourOfPair> sortedLatticeSitesEncoded;
  sortedLatticeSitesEncoded.reserve(numSortedSites);

  unordered_map<size_t, size_t> latticeIdToSortedIndexMap;
  latticeIdToSortedIndexMap.reserve(numSortedSites);

  size_t index = 0;

  // Iterate over all the sites and get store the sorted lattice Ids
  for (size_t site1 = 0; site1 < numLattices; site1++)
  {
    for (const size_t site2 : smallConfig.GetNeighborLatticeIdVectorOfLattice(site1, 1))
    {
      pair<size_t, size_t> canonicalLatticeIdJumpPair = (site1 < site2)
                                                            ? make_pair(site1, site2)
                                                            : make_pair(site2, site1);

      if (pairToIndexHashMap_.find(canonicalLatticeIdJumpPair) != pairToIndexHashMap_.end())
      {
        // Pair already exists
        continue;
      }

      // Encode a jump Pair to a index of a vector
      pairToIndexHashMap_[canonicalLatticeIdJumpPair] = index;
      index += 1;

      auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForPair(
          smallConfig,
          canonicalLatticeIdJumpPair,
          maxBondOrder_,
          canonicalReferenceMap_,
          symmetryOperations_);

      // Now create the latticeId
      latticeIdToSortedIndexMap.clear();
      for (size_t i = 0; i < numSortedSites; i++)
      {
        latticeIdToSortedIndexMap[canonicalSortedLatticeIds[i]] = i;
      }

      // now get the neighbourlist for the latticeIdPair from the tiledSupercell
      // and build the canonicalSortedLatticeAndConfigIdVector
      // <nnLatticeId, cubeIdx>
      sortedLatticeSitesEncoded.clear();
      sortedLatticeSitesEncoded = tiledSupercell.GetNeighboringLatticeIdSetOfPair(
          canonicalLatticeIdJumpPair,
          maxBondOrder_,
          removeLatticeIdDuplicates);

      std::sort(sortedLatticeSitesEncoded.begin(), sortedLatticeSitesEncoded.end(),
                [&](const auto &a, const auto &b)
                {
                  return latticeIdToSortedIndexMap.at(
                             a.latticeSiteInfo.latticeId) < latticeIdToSortedIndexMap.at(b.latticeSiteInfo.latticeId);
                });

      symmetricLatticeSiteEncodedMappings_.emplace_back(move(sortedLatticeSitesEncoded));
    }
  }
}

void KRAPredictorTLMC::PrintKRAPredictorInfo() const
{
  const int width = 80;
  const int labelWidth = 40;
  const int valueWidth = width - labelWidth - 2; // 2 for spacing

  // Header
  cout << string(width, '-') << "\n";
  cout << setw((width + 22) / 2) << right << "KRA Predictor Info" << "\n";
  cout << string(width, '-') << "\n";

  // Table of values
  cout << left << setw(labelWidth) << "Max bond order:"
       << right << setw(valueWidth) << maxBondOrder_ << "\n";

  cout << left << setw(labelWidth) << "Max cluster size:"
       << right << setw(valueWidth) << maxClusterSize_ << "\n";

  cout << left << setw(labelWidth) << "Reference jump direction:"
       << right << setw(valueWidth - 3)
       << "(" << referenceJumpDirection_.transpose() << ")" << "\n";

  cout << left << setw(labelWidth) << "Canonical reference map size:"
       << right << setw(valueWidth) << canonicalReferenceMap_.size() << "\n";

  cout << left << setw(labelWidth) << "Number of symmetry operations:"
       << right << setw(valueWidth) << symmetryOperations_.size() << "\n";

  cout << left << setw(labelWidth) << "Encoded orbits for pair size:"
       << right << setw(valueWidth) << encodedOrbitsForPair_.size() << "\n";

  cout << left << setw(labelWidth) << "Size of precomputed vector map for jump pairs:"
       << right << setw(valueWidth - 10) << "(" << pairToIndexHashMap_.size() << ", "
       << symmetricLatticeSiteEncodedMappings_.size() << ")" << "\n";

  cout << left << setw(labelWidth) << "Kinetic ECIs map size:"
       << right << setw(valueWidth) << KECIsMap_.size() << "\n";

  cout << left << setw(labelWidth) << "KECIs per element:" << "\n";
  for (const auto &[element, vec] : KECIsMap_)
  {
    cout << "  " << left << setw(labelWidth - 2)
         << element.GetElementString()
         << right << setw(valueWidth)
         << vec.size() << "\n";
  }

  cout << string(width, '-') << "\n\n";
}