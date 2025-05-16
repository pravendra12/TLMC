/*******************************************************************************************
 * Copyright (c) 2025. All rights reserved.                                     *
 * @Author: Pravendra Patel                                                                *
 * @Date:    2025-05-08                              *
 * @Last Modified by: Pravendra Patel                                                      *
 * @Last Modified time: 2025-05-08                   *
 *******************************************************************************************/

#include "KRAPredictor.h"
#include <chrono>

using namespace chrono;

#include <iostream>
#include <map>
#include <vector>

size_t estimateVectorSize(const std::vector<std::vector<size_t>>& v) {
    size_t total = sizeof(v);  // vector of vectors
    total += v.capacity() * sizeof(std::vector<size_t>); // inner vectors metadata

    for (const auto& inner : v) {
        total += sizeof(inner);  // each inner vector
        total += inner.capacity() * sizeof(size_t);  // actual data
    }

    return total;
}

size_t estimateMapSize(const std::map<string, std::vector<std::vector<size_t>>>& orbitMap) {
    size_t total = 0;
    size_t nodeOverhead = 48;  // map node overhead (approximate)

    for (const auto& [key, value] : orbitMap) {
        total += nodeOverhead;
        total += sizeof(key);   // key size
        total += estimateVectorSize(value); // value size
    }

    return total;
}


KRAPredictor::KRAPredictor(
    const string &predictorFilename,
    const Config &config,
    const set<Element> &elementSet) : betaKRA_W_(ReadParametersFromJson(predictorFilename,
                                                                     "EKRA_W", "beta_W")),
                                   interceptKRA_W_(ReadParametersFromJson(
                                       predictorFilename,
                                       "EKRA_W", "intercept_W")(0)),
                                   betaKRA_Ta_(ReadParametersFromJson(
                                       predictorFilename,
                                       "EKRA_Ta", "beta_Ta")),
                                   interceptKRA_Ta_(ReadParametersFromJson(
                                       predictorFilename,
                                       "EKRA_Ta", "intercept_Ta")(0)),
                                   elementSet_([&]()
                                               {
                                                set<Element> cleanedSet = elementSet;
                                                cleanedSet.erase(Element("X"));  
                                                return cleanedSet; }()),
                                   maxBondOrder_(
                                       ReadParameterFromJson(
                                           predictorFilename,
                                           "maxBondOrderKRA")),
                                   maxBondOrderOfCluster_(
                                       ReadParameterFromJson(
                                           predictorFilename,
                                           "maxBondOrderOfClusterKRA")),
                                   maxClusterSize_(
                                       ReadParameterFromJson(
                                           predictorFilename,
                                           "maxClusterSizeKRA")),
                                   equivalentSiteEncoding_(
                                       GetEquivalentSiteEncoding3BarSymmetry(
                                           config,
                                           maxBondOrder_)),
                                   latticePairToSSVectorMap_(
                                       GetSymmetricallySortedLatticePairMap(
                                           config,
                                           maxBondOrder_))

{
  // Check for maxBondOrder, maxClusterSize and maxBondOrderOfCluster using the
  // json file

  cout << "Max Bond Order for KRA: " << maxBondOrder_ << endl;
  cout << "Max Bond Order of Clusters for KRA: " << maxBondOrderOfCluster_ << endl;
  cout << "Max Cluster Size for KRA:  " << maxClusterSize_ << endl;

}

double KRAPredictor::GetKRA(const Config &config,
                            const pair<size_t, size_t> &latticeIdJumpPair) const
{

  const auto &elementFirst = config.GetElementOfLattice(latticeIdJumpPair.first);
  const auto &elementSecond = config.GetElementOfLattice(latticeIdJumpPair.second);

  string migratingAtom = (elementFirst.GetElementString() == "X")
                             ? elementSecond.GetElementString()
                             : elementFirst.GetElementString();

  // Canonical Pair
  size_t id1 = latticeIdJumpPair.first;
  size_t id2 = latticeIdJumpPair.second;

  pair<size_t, size_t> canonicalJumpPair = (id1 < id2)
                                               ? std::make_pair(id1, id2)
                                               : std::make_pair(id2, id1);

  // O(1)
  auto ssVector = latticePairToSSVectorMap_.at(canonicalJumpPair);

  // Orbit Map
  auto startOrbitMap = high_resolution_clock::now();

  auto orbitMap = GetOrbits(config,
                            maxClusterSize_,
                            maxBondOrderOfCluster_,
                            ssVector,
                            equivalentSiteEncoding_);

  cout << "Size of a single map: " << estimateMapSize(orbitMap) / (1024.0 * 1024.0) << " MB\n" << endl;

  auto endOrbitMap = high_resolution_clock::now();

  auto durationOrbitMap = duration_cast<microseconds>(endOrbitMap - startOrbitMap);
  // cout << "Time to compute orbit map: " << durationOrbitMap.count() << " microseconds" << endl;

  // string basisType = "Occupation";

  // Get the LCE
  auto startLCE = high_resolution_clock::now();

  VectorXd localEnvEncoding = GetLocalEnvironmentEncoding(config,
                                                          elementSet_,
                                                          basisType_,
                                                          orbitMap)
                                  .transpose();

  auto endLCE = high_resolution_clock::now();

  auto durationLCE = duration_cast<microseconds>(endLCE - startLCE);
  // cout << "Time to compute LCE: " << durationLCE.count() << " microseconds" << endl;

  if (localEnvEncoding.size() != betaKRA_W_.size())
  {
    cerr << "Error in GetKRA: Mismatch in local environment encoding size!" << endl;
    cerr << "Expected size: " << betaKRA_W_.size() << " (from betaKRA_W_)" << endl;
    cerr << "Actual size  : " << localEnvEncoding.size() << " (from localEnvEncoding)" << endl;

    cerr << "This likely indicates an issue with how the local environment was encoded." << endl;
    cerr << "Check if all expected elements are present and properly mapped." << endl;

    cerr << "Details for debugging:" << endl;
    cerr << " - Size of betaKRA_Ta_: " << betaKRA_Ta_.size() << endl;
    cerr << " - Size of betaKRA_W_ : " << betaKRA_W_.size() << endl;
    cerr << " - Element Set: ";
    for (const auto &ele : elementSet_)
    {
      cerr << ele.GetElementString() << " ";
    }
    cerr << endl;

    exit(2);
  }

  double eKRA;

  if (migratingAtom == "W")
  {
    eKRA = betaKRA_W_.dot(localEnvEncoding) + interceptKRA_W_;
  }
  else if (migratingAtom == "Ta")
  {
    eKRA = betaKRA_Ta_.dot(localEnvEncoding) + interceptKRA_Ta_;
  }
  else
  {
    cout << "Fitting coefficients not found for " << migratingAtom << endl;
    exit(1);
  }

  // cout << "Coeff Ta : " << betaKRA_Ta_ << endl;
  // cout << "Intercept Ta: " << interceptKRA_Ta_ << endl;
  // cout << "Coeff W : " << betaKRA_W_ << endl;

  // cout << "EKRA : " << eKRA << endl;

  return eKRA;
}

// Precompute all the symmetrically sorted pairs
static PairMap GetSymmetricallySortedLatticePairMap(
    const Config &config,
    const size_t maxBondOrder)
{
  PairMap ssVectorLatticePairMap;

  auto numLattice = config.GetNumLattices();

  for (size_t id1 = 0; id1 < numLattice; id1++)
  {
    auto id1FirstNN = config.GetNeighborLatticeIdVectorOfLattice(id1, 1);

    for (auto &id2 : id1FirstNN)
    {
      // Canonical Pair
      std::pair<size_t, size_t> canonicalJumpPair = (id1 < id2)
                                                        ? std::make_pair(id1, id2)
                                                        : std::make_pair(id2, id1);

      // Check if lattice pairs are already present in the global map (optional, depending on your need)
      if (ssVectorLatticePairMap.find(canonicalJumpPair) == ssVectorLatticePairMap.end())
      {
        // Get symmetrically sorted vector for the lattice pair under 3 Bar symmetry
        auto ssVector = GetSSVector3FSymmetry(
            config,
            canonicalJumpPair,
            maxBondOrder);

        ssVectorLatticePairMap[canonicalJumpPair] = ssVector;
      }
    }
  }

  return ssVectorLatticePairMap;
}