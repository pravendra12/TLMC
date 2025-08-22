/*******************************************************************************
 * Copyright (c) 2025. All rights reserved.
 * @Author: Pravendra Patel
 * @Date:    2025-06-01
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 *******************************************************************************/

/*! @file KRAPredictor.cpp
    @brief File contains implementation of Kinetically Resolved Activation Barrier
 */

#include "KRAPredictor.h"

KRAPredictor::KRAPredictor(
    const string &predictorFilename,
    const Config &referenceConfig,
    const set<Element> &elementSet) : betaKRA_W_(ReadParametersFromJson(predictorFilename,
                                                                        "EKRA_W", "beta_W")),
                                      interceptKRA_W_(
                                          ReadParametersFromJson(
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
                                              referenceConfig,
                                              maxBondOrder_)),
                                      latticePairToSSVectorMap_(
                                          GetSymmetricallySortedLatticePairMap(
                                              referenceConfig,
                                              maxBondOrder_))

{

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

  auto orbitMap = GetOrbits(config,
                            maxClusterSize_,
                            maxBondOrderOfCluster_,
                            ssVector,
                            equivalentSiteEncoding_);

  // Get the LCE
  VectorXd localEnvEncoding = GetLocalEnvironmentEncoding(
                                  config,
                                  elementSet_,
                                  basisType_,
                                  orbitMap)
                                  .transpose();

  // cout << localEnvEncoding.size() << endl;

  if (localEnvEncoding.size() != betaKRA_W_.size())
  {
    std::ostringstream oss;
    oss << "Error in KRAPredictor: Mismatch in local environment encoding size!\n"
        << "Expected size: " << betaKRA_W_.size() << " (from betaKRA_W_)\n"
        << "Actual size  : " << localEnvEncoding.size() << " (from localEnvEncoding)\n\n"
        << "This likely indicates an issue with how the local environment was encoded.\n"
        << "Check if all expected elements are present and properly mapped.\n\n"
        << "Details for debugging:\n"
        << " - Size of betaKRA_Ta_: " << betaKRA_Ta_.size() << "\n"
        << " - Size of betaKRA_W_ : " << betaKRA_W_.size() << "\n"
        << " - Element Set: ";

    for (const auto &ele : elementSet_)
    {
      oss << ele.GetElementString() << " ";
    }

    throw std::invalid_argument(oss.str());
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
    throw std::invalid_argument("Error in KRAPredictor: Fitting coefficients not found for atom: " + migratingAtom);
  }

  return eKRA;
}

PairMap GetSymmetricallySortedLatticePairMap(
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