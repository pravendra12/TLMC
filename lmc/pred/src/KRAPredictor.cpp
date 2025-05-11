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

// Precompute all the symmetrically sorted pairs

/*
using PairMap = unordered_map<pair<size_t, size_t>,
                              vector<size_t>,
                              boost::hash<pair<size_t, size_t>>>;

PairMap GetSymmetricallySortedLatticePairMap(const Config &config,
                                             const size_t maxBondOrder)
{
  PairMap ssVectorLatticePairMap;

  auto numLattice = config.GetNumLattices();

  for (size_t id1 = 0; id1 < numLattice; id1++)
  {
    auto id1FirstNN = config.GetNeighborLatticeIdVectorOfLattice(id1, 1);

    for (auto &id2 : id1FirstNN)
    {
      pair<size_t, size_t> latticePair1 = {id1, id2};
      pair<size_t, size_t> latticePair2 = {id2, id1};

      // Check if lattice pairs are already present in the global map (optional, depending on your need)
      if (ssVectorLatticePairMap.find(latticePair1) == ssVectorLatticePairMap.end() &&
          ssVectorLatticePairMap.find(latticePair2) == ssVectorLatticePairMap.end())
      {
        // Get symmetrically sorted vector for the lattice pairs
        auto ssVector1 = GetSortedLatticeStatesForPairUnder3BarSymmetry(
            config,
            latticePair1,
            maxBondOrder);

        auto ssVector2 = GetSortedLatticeStatesForPairUnder3BarSymmetry(
            config,
            latticePair2,
            maxBondOrder);

        ssVectorLatticePairMap[latticePair1] = ssVector1;
        ssVectorLatticePairMap[latticePair2] = ssVector2;
      }
    }
  }

  return ssVectorLatticePairMap;
}
*/

KRAPredictor::KRAPredictor(
    const string &predictorFilename,
    const Config &config,
    const set<Element> &elementSet,
    const size_t maxBondOrder,
    const size_t maxClusterSize) : betaKRA_W_(ReadParametersFromJson(predictorFilename,
                                                                     "EKRA_W", "beta_W")),
                                   interceptKRA_W_(ReadParametersFromJson(predictorFilename,
                                                                          "EKRA_W", "intercept_W")(0)),
                                   betaKRA_Ta_(ReadParametersFromJson(predictorFilename,
                                                                      "EKRA_Ta", "beta_Ta")),
                                   interceptKRA_Ta_(ReadParametersFromJson(predictorFilename,
                                                                           "EKRA_Ta", "intercept_Ta")(0)),
                                   elementSet_([&]()
                                               {
                                                set<Element> cleanedSet = elementSet;
                                                cleanedSet.erase(Element("X"));  
                                                return cleanedSet; }()),
                                   maxBondOrder_(maxBondOrder),
                                   maxClusterSize_(maxClusterSize)
{
}

double KRAPredictor::GetKRA(const Config &config,
                            const pair<size_t, size_t> &latticeIdJumpPair) const
{

  const auto &elementFirst = config.GetElementOfLattice(latticeIdJumpPair.first);
  const auto &elementSecond = config.GetElementOfLattice(latticeIdJumpPair.second);

  string migratingAtom = (elementFirst.GetElementString() == "X")
                             ? elementSecond.GetElementString()
                             : elementFirst.GetElementString();

  // cout << "Migrating Atom : " << migratingAtom << endl;

  // Equivalent encoding
  // Need to update the below funciton as in principle it should be independent of
  // latticeIdJumpPair
  auto startEqSites3Bar = high_resolution_clock::now();

  auto eqSites3Bar = GetEquivalentSitesUnder3BarSymmetry(config,
                                                         latticeIdJumpPair,
                                                         maxBondOrder_);
  // print2DVector(equivalentSites3Bar_);
  auto endEqSites3Bar = high_resolution_clock::now();

  auto durationEqSites = duration_cast<microseconds>(endEqSites3Bar - startEqSites3Bar);
  // cout << "Time to compute eq sites: " << durationEqSites.count() << " microseconds" << endl;

  // print2DVector(eqSites3Bar);

  // OrbitMap
  // size_t maxClusterSize = 3;
  auto startOrbitMap = high_resolution_clock::now();
  auto orbitMap = GetOrbits(config,
                            maxClusterSize_,
                            maxBondOrder_,
                            eqSites3Bar);

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


  // cout << "Size of localEnvEncoding : " << localEnvEncoding.size() << endl;
  // cout << "Size of beta Ta : " << betaKRA_Ta_.size() << endl;
  // cout << "Size of beta W : " << betaKRA_W_.size() << endl;

  if (localEnvEncoding.size() != 414)
  {
    cout << "The size of localEnvEncoding is not 414, issue in GetKRA!" << endl;
    for (auto ele : elementSet_)
    {
      cout << ele.GetElementString() << endl;
    }
    exit(2);
  }

  double eKRA;

  // cout << betaKRA_Ta_.size() << endl;
  // cout << localEnvEncoding.size() << endl;

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
