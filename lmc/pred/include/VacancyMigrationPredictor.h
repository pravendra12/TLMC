#ifndef LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_
#define LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_

#include "VacancyMigrationBarrierPredictor.h"
#include "PotentialEnergyEstimator.h"

using namespace std;

class VacancyMigrationPredictor
{
public:
  VacancyMigrationPredictor(
      const string &predictorFilename,
      const Config &referenceConfig,
      const Config &supercellConfig,
      const set<Element> &elementSet,
      size_t maxClusterSize,
      size_t maxBondOrder);

  // In current model 
  // First lattice Id in the latticeIdJumpPair is the migrating atom Id based on
  // which the symmetrically sorted vector is computed and it is pivoted towards
  // the migrating atom Id.
  
  // latticeIdJumpPair : <NeighbourAtomId, VacancyId>
  pair<double, double> getBarrierAndEnergyChange(
      const Config &config,
      const pair<size_t, size_t> &latticeIdJumpPair) const;

private:
  VacancyMigrationBarrierPredictor barrierPredictor_;
  PotentialEnergyEstimator energyChangePredictor_;
};

#endif // LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_
