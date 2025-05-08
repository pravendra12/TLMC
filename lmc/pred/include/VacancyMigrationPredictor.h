#ifndef LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_
#define LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_

#include "KRAPredictor.h"
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

  pair<double, double> GetBarrierAndDeltaE(
      const Config &config,
      const pair<size_t, size_t> &latticeIdJumpPair) const;

private:
  KRAPredictor eKRAPredictor;
  PotentialEnergyEstimator energyChangePredictor_;
};

#endif // LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_
