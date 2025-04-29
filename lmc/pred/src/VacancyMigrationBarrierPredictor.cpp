#include "VacancyMigrationBarrierPredictor.h"

// This version of the barrier predictor uses the upto 3rd nearest neighbours
// and also uses 6 Fold symmetry encoding, for more details see the docs ..
// Lattice Jump pair = {vacancyId, migratingAtomId} for forward barrier predicton

VacancyMigrationBarrierPredictor::VacancyMigrationBarrierPredictor(
    const Config &config,
    const set<Element> &elementSet,
    const string &predictorFilename) : betaBarrier_(ReadParametersFromJson(predictorFilename,
                                                                           "barrier",
                                                                           "beta_barrier")),
                                       meanXEncoding_(ReadParametersFromJson(predictorFilename,
                                                                             "barrier",
                                                                             "mean_X_features")),
                                       stdXEncoding_(ReadParametersFromJson(predictorFilename,
                                                                            "barrier",
                                                                            "std_X_features")),
                                       meanYBarrier_(ReadParametersFromJson(predictorFilename,
                                                                            "barrier",
                                                                            "mean_Y_barrier")(0)),
                                       stdYBarrier_(ReadParametersFromJson(predictorFilename,
                                                                           "barrier",
                                                                           "std_Y_barrier")(0)),

                                       // barrier_fitted_parameters_(ReadParametersFromJson(predictorFilename,
                                       //                                                   "barrier")),
                                       // adjusted_beta_barrier_(barrier_fitted_parameters_.first),
                                       // adjusted_intercept_barrier_(barrier_fitted_parameters_.second),
                                       oneHotEncodingMap_(
                                           GetOneHotEncodeHashmap(elementSet)),
                                       encodingKFoldRotation_(
                                           GetEquivalentSitesUnderKFoldRotation(config,
                                                                                maxBondOrderForBarrierPrediction_,
                                                                                kFoldRotation_)),
                                       symmetricallySortedVectorMap_(
                                           ComputeSymmetricallySortedVectorMap(config,
                                                                               maxBondOrderForBarrierPrediction_))
{
}

/*
double VacancyMigrationBarrierPredictor::GetBarrier(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair) const
{
  Element migratingAtom;
  // Just a check
  if (config.GetElementOfLattice(latticeIdJumpPair.first) == Element("X"))
  {
    migratingAtom = config.GetElementOfLattice(latticeIdJumpPair.second);
  }
  else
  {
    migratingAtom = config.GetElementOfLattice(latticeIdJumpPair.first);
  }

  // pair<size_t, size_t> pair_test = {0, 1};
  // O(1)
  auto sortedLatticeVector = symmetricallySortedVectorMap_.at(latticeIdJumpPair);

  // auto startGetEncodingVector = std::chrono::high_resolution_clock::now();

  // encoding for given lattice pair and the migrating atom
  VectorXd migratingAtomEncodingVector = GetEncodingMigratingAtomPair(config,
                                                                      encoding3FoldRotation_,
                                                                      sortedLatticeVector,
                                                                      oneHotEncodingMap_,
                                                                      migratingAtom);

  // cout << migratingAtomEncodingVector.transpose() << endl;

  // auto endGetEncodingVector = std::chrono::high_resolution_clock::now();

  // std::chrono::duration<double> duration = endGetEncodingVector - startGetEncodingVector;

  // cout << "Time taken by GetEncodingMigrationAtomPair computation after changes: " << duration.count() << endl;

  double barrier = migratingAtomEncodingVector.dot(adjusted_beta_barrier_) +
                   adjusted_intercept_barrier_;

  return barrier;
}
*/

double VacancyMigrationBarrierPredictor::GetBarrier(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair) const
{
  // latticeIdJumpPair : {vacancyId, migratingAtomId}

  Element migratingAtom;

  // Just a check
  if (config.GetElementOfLattice(latticeIdJumpPair.first) == Element("X"))
  {
    migratingAtom = config.GetElementOfLattice(latticeIdJumpPair.second);
  }
  else if (config.GetElementOfLattice(latticeIdJumpPair.second) == Element("X"))
  {
    migratingAtom = config.GetElementOfLattice(latticeIdJumpPair.first);
  }
  else
  {
    std::cout << "No vacancy found in the given Jump Pair" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // pair<size_t, size_t> pair_test = {0, 1};
  // O(1)
  auto sortedLatticeVector = symmetricallySortedVectorMap_.at(latticeIdJumpPair);

  // auto startGetEncodingVector = std::chrono::high_resolution_clock::now();

  // encoding for given lattice pair and the migrating atom
  VectorXd barrierEncoding = GetEncodingMigratingAtomPair(config,
                                                          encodingKFoldRotation_,
                                                          sortedLatticeVector,
                                                          oneHotEncodingMap_,
                                                          migratingAtom);

  // cout << migratingAtomEncodingVector.transpose() << endl;

  // auto endGetEncodingVector = std::chrono::high_resolution_clock::now();

  // std::chrono::duration<double> duration = endGetEncodingVector - startGetEncodingVector;

  // cout << "Time taken by GetEncodingMigrationAtomPair computation after changes: " << duration.count() << endl;

  // Scale the feature vector

  VectorXd scaledEncodings = (barrierEncoding.array() - meanXEncoding_.array()) / stdXEncoding_.array();

  double scaledBarrier = scaledEncodings.dot(betaBarrier_);

  // Inverse scaling to get the correct barrier

  double barrier = scaledBarrier * stdYBarrier_ + meanYBarrier_;

  return barrier;
}

PairMap ComputeSymmetricallySortedVectorMap(const Config &config,
                                            const size_t maxBondOrder)
{

  PairMap symmetricallySortedVectorMap;
  size_t numLattices = config.GetNumLattices();

// Parallel loop with thread-local storage for the local map
#pragma omp parallel
  {
    // Thread-local map to store results for each thread
    unordered_map<pair<size_t, size_t>, vector<size_t>, boost::hash<pair<size_t, size_t>>> localMap;

#pragma omp for
    for (size_t id1 = 0; id1 < numLattices; ++id1)
    {
      // First nearest neighbors of id1
      auto id1FirstNN = config.GetNeighborLatticeIdVectorOfLattice(id1, 1);

      for (auto &id2 : id1FirstNN)
      {
        pair<size_t, size_t> latticePair1 = {id1, id2};
        pair<size_t, size_t> latticePair2 = {id2, id1};

        // Check if lattice pairs are already present in the global map (optional, depending on your need)
        if (symmetricallySortedVectorMap.find(latticePair1) == symmetricallySortedVectorMap.end() &&
            symmetricallySortedVectorMap.find(latticePair2) == symmetricallySortedVectorMap.end())
        {
          // Get symmetrically sorted vector for the lattice pairs
          auto symmetricallySortedVector1 = GetSortedLatticeVectorStateOfPair(config, latticePair1, maxBondOrder);
          auto symmetricallySortedVector2 = GetSortedLatticeVectorStateOfPair(config, latticePair2, maxBondOrder);

          // Add both directions to the thread-local map
          localMap[latticePair1] = symmetricallySortedVector1;
          localMap[latticePair2] = symmetricallySortedVector2;
        }
      }
    }

// Lock and merge thread-local map into global map after the parallel block
#pragma omp critical
    {
      symmetricallySortedVectorMap.insert(localMap.begin(), localMap.end());
    }
  }

  return symmetricallySortedVectorMap;
}
