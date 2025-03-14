#include "VacancyMigrationBarrierPredictor.h"

#include "PrintUtility.h"


VacancyMigrationBarrierPredictor::VacancyMigrationBarrierPredictor(
    const Config &config, 
    const set<Element> &elementSet, 
    const size_t &maxBondOrder,
    const string &predictorFilename): 
                barrier_fitted_parameters_(
                     ReadParametersFromJson(predictorFilename, 
                                            "barrier")),
                adjusted_beta_barrier_(barrier_fitted_parameters_.first),
                adjusted_intercept_barrier_(barrier_fitted_parameters_.second),
                oneHotEncodingMap_(
                  GetOneHotEncodeHashmap(elementSet)),
                encoding3FoldRotation_(
                  GetEquivalentSites3Fold(config, 
                                          maxBondOrderForBarrierPrediction_))                
                // maxBondOrderForBarrierPrediction_(maxBondOrder)
{
  ComputeSymmetricallySortedVectorMap(config, 
                                      maxBondOrderForBarrierPrediction_, 
                                      symmetricallySortedVectorMap_);
}

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

  auto sortedLatticeVector = symmetricallySortedVectorMap_.at(latticeIdJumpPair);

  // print1DVector(sortedLatticeVector);
  
  // encoding for given lattice pair and the migrating atom
  VectorXd migratingAtomEncodingVector = GetEncodingMigratingAtomPair(config, 
                                                             encoding3FoldRotation_, 
                                                             sortedLatticeVector, 
                                                             oneHotEncodingMap_, 
                                                             migratingAtom);
  // print2DVector(encoding3FoldRotation_);
  // cout << migratingAtomEncodingVector.transpose() << endl;


  double barrier = migratingAtomEncodingVector.dot(adjusted_beta_barrier_) + 
                   adjusted_intercept_barrier_;



  return barrier;
}


void ComputeSymmetricallySortedVectorMap(const Config &config, 
  const size_t maxBondOrder, 
  unordered_map<pair<size_t, size_t>, vector<size_t>, 
  boost::hash<pair<size_t, size_t>>> &symmetricallySortedVectorMap) 
{
  size_t numLattices = config.GetNumLattices();
  mutex mapMutex; // Mutex for thread-safe map insertion

  #pragma omp parallel for
  for (size_t id1 = 0; id1 < numLattices; ++id1) 
  {
    // First nearest neighbors of id1
    auto id1FirstNN = config.GetNeighborLatticeIdVectorOfLattice(id1, 1);
    
    // local map
    unordered_map<pair<size_t, size_t>, 
                  vector<size_t>, 
                  boost::hash<pair<size_t, size_t>>> localMap;

    for (auto &id2 : id1FirstNN) 
    {
      pair<size_t, size_t> latticePair = {id1, id2};

      auto symmetricallySortedVector = GetSortedLatticeVectorStateOfPair(config, 
                                                                         latticePair, 
                                                                         maxBondOrder);
      // Appending to the local map
      localMap[latticePair] = symmetricallySortedVector;
    }

    // Lock and merge local map into global map
    lock_guard<mutex> lock(mapMutex);
    symmetricallySortedVectorMap.insert(localMap.begin(), localMap.end());
  }
}
