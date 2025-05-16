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

size_t estimateVectorSize(const std::vector<std::vector<size_t>> &v)
{
  size_t total = sizeof(v);                            // vector of vectors
  total += v.capacity() * sizeof(std::vector<size_t>); // inner vectors metadata

  for (const auto &inner : v)
  {
    total += sizeof(inner);                     // each inner vector
    total += inner.capacity() * sizeof(size_t); // actual data
  }

  return total;
}

template <typename T>
size_t estimateMapSize(const std::map<T, std::vector<std::vector<size_t>>> &orbitMap)
{
  size_t total = 0;
  size_t nodeOverhead = 48; // map node overhead (approximate)

  for (const auto &[key, value] : orbitMap)
  {
    total += nodeOverhead;
    total += sizeof(key);               // key size
    total += estimateVectorSize(value); // value size
  }

  return total;
}

template <typename T>
size_t estimateVectorPairSize(const std::vector<std::pair<std::vector<T>, std::vector<std::vector<T>>>> &orbitVec)
{
  size_t total = 0;

  for (const auto &[keyVec, valueVecVec] : orbitVec)
  {
    // Approximate overhead per vector object
    size_t vectorOverhead = 24; // typical small vector metadata size (capacity, size, pointer)

    // Size of key vector content
    total += vectorOverhead + (keyVec.capacity() * sizeof(T));

    // Size of value vector of vectors content
    total += vectorOverhead; // outer vector overhead
    for (const auto &innerVec : valueVecVec)
    {
      total += vectorOverhead + (innerVec.capacity() * sizeof(T));
    }
  }

  return total;
}

#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <utility>
#include <boost/functional/hash.hpp>

// Estimate memory usage of a vector
template <typename T>
size_t estimateVectorSize(const std::vector<T> &vec)
{
  return sizeof(vec) + sizeof(T) * vec.capacity();
}

// Estimate memory usage of a nested vector
template <typename T>
size_t estimateNestedVectorSize(const std::vector<std::vector<T>> &vec)
{
  size_t total = sizeof(vec);
  for (const auto &inner : vec)
  {
    total += estimateVectorSize(inner);
  }
  return total;
}

// Estimate memory usage of map<string, vector<vector<size_t>>>
size_t estimateInnerMapSize(const std::map<std::string, std::vector<std::vector<size_t>>> &innerMap)
{
  size_t total = sizeof(innerMap);
  for (const auto &[key, val] : innerMap)
  {
    total += sizeof(key) + key.capacity();  // key string
    total += estimateNestedVectorSize(val); // value vector<vector<size_t>>
  }
  return total;
}

// Estimate memory usage of PairToOrbitMap
size_t estimatePairToOrbitMapSize(const std::unordered_map<std::pair<size_t, size_t>,
                                                           std::map<std::string, std::vector<std::vector<size_t>>>,
                                                           boost::hash<std::pair<size_t, size_t>>> &pairToOrbitMap)
{
  size_t total = sizeof(pairToOrbitMap);
  size_t nodeOverhead = 32; // Approximate overhead per unordered_map entry

  for (const auto &[key, innerMap] : pairToOrbitMap)
  {
    total += nodeOverhead;
    total += sizeof(key);                    // the pair key
    total += estimateInnerMapSize(innerMap); // the value map<string, vector<vector<size_t>>>
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
                                      latticePairToOrbitMap_(
                                          GetPairToOrbitMap(config))

{
  // Check for maxBondOrder, maxClusterSize and maxBondOrderOfCluster using the
  // json file

  cout << "Max Bond Order for KRA: " << maxBondOrder_ << endl;
  cout << "Max Bond Order of Clusters for KRA: " << maxBondOrderOfCluster_ << endl;
  cout << "Max Cluster Size for KRA:  " << maxClusterSize_ << endl;

  size_t bytes = estimatePairToOrbitMapSize(latticePairToOrbitMap_);
  std::cout << "Memory usage: " << bytes / (1024.0 * 1024.0) << " MB\n";
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
  // auto ssVector = latticePairToSSVectorMap_.at(canonicalJumpPair);

  // Orbit Map
  //  auto startOrbitMap = high_resolution_clock::now();
  //
  //  auto orbitMap = GetOrbits(config,
  //                            maxClusterSize_,
  //                            maxBondOrderOfCluster_,
  //                            ssVector,
  //                            equivalentSiteEncoding_);
  //
  //  cout << "Size of a single map: " << estimateMapSize(orbitMap) / (1024.0 * 1024.0) << " MB\n"
  //       << endl;
  //
  //  auto endOrbitMap = high_resolution_clock::now();
  //
  //  auto durationOrbitMap = duration_cast<microseconds>(endOrbitMap - startOrbitMap);
  // cout << "Time to compute orbit map: " << durationOrbitMap.count() << " microseconds" << endl;

  // string basisType = "Occupation";

  auto startOrbitMap = high_resolution_clock::now();
  auto orbitMap = latticePairToOrbitMap_.at(canonicalJumpPair);
  auto endOrbitMap = high_resolution_clock::now();
  auto durationOrbitMap = duration_cast<microseconds>(endOrbitMap - startOrbitMap);

  cout << "Time to compute orbit map: " << durationOrbitMap.count() << " microseconds" << endl;

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

PairToOrbitMap KRAPredictor::GetPairToOrbitMap(const Config &config)
{
  PairToOrbitMap latticePairToOrbitMap;

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
      if (latticePairToOrbitMap.find(canonicalJumpPair) == latticePairToOrbitMap.end())
      {
        // Get symmetrically sorted vector for the lattice pair under 3 Bar symmetry
        auto ssVector = GetSSVector3FSymmetry(
            config,
            canonicalJumpPair,
            maxBondOrder_);

        auto orbitMap = GetOrbits(config,
                                  maxClusterSize_,
                                  maxBondOrderOfCluster_,
                                  ssVector,
                                  equivalentSiteEncoding_);

        latticePairToOrbitMap[canonicalJumpPair] = orbitMap;
      }
    }
  }

  return latticePairToOrbitMap;
}
