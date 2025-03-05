#include "VacancyMigrationBarrierPredictor.h"

VacancyMigrationBarrierPredictor::VacancyMigrationBarrierPredictor(
    const Config &config, 
    const std::set<Element> &elementSet, 
    const size_t &maxBondOrder, 
    const std::string &predictorFilename): 
                barrier_fitted_parameters_(
                     ReadParametersFromJson(predictorFilename, 
                                            "barrier")),
                adjusted_beta_barrier_(barrier_fitted_parameters_.first),
                adjusted_intercept_barrier_(barrier_fitted_parameters_.second),
                oneHotEncodingMap_(GetOneHotEncodeHashmap(move(elementSet))),
                encoding3FoldRotation_(GetEquivalentSites3Fold(move(config), 
                                                               maxBondOrder)),
                // For now only works for bond order of 2
                maxBondOrder_(2)
{
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

  auto sortedLatticeVector = GetSortedLatticeVectorStateOfPair(config, 
                                                               latticeIdJumpPair, 
                                                               maxBondOrder_);
  // encoding for given lattice pair and the migrating atom
  VectorXd migratingAtomEncodingVector = GetEncodingMigratingAtomPairs(config, 
                                                             encoding3FoldRotation_, 
                                                             sortedLatticeVector, 
                                                             oneHotEncodingMap_, 
                                                             migratingAtom);


  double barrier = migratingAtomEncodingVector.dot(adjusted_beta_barrier_) + 
                   adjusted_intercept_barrier_;



  return barrier;
}
