#include "VacancyMigrationPredictor.h"

VacancyMigrationPredictor::VacancyMigrationPredictor(
    const string &predictorFilename,
    const Config &referenceConfig,
    const Config &supercellConfig,
    const set<Element> &elementSet,
    size_t maxClusterSize,
    size_t maxBondOrder) : eKRAPredictor(predictorFilename,
                                         referenceConfig,
                                         elementSet),
                           energyChangePredictor_(predictorFilename,
                                                  referenceConfig,
                                                  supercellConfig,
                                                  elementSet)
{
}

// (barrier, dE)
pair<double, double> VacancyMigrationPredictor::GetBarrierAndDeltaE(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair) const
{
  double eKRA = eKRAPredictor.GetKRA(config, 
                                     latticeIdJumpPair);

  double dE = energyChangePredictor_.GetDeThreadSafe(config, 
                                                     latticeIdJumpPair);
  
  // EKRA = Ea - 1/2*dE

  double barrier = eKRA + (dE/2);
  
  return pair<double, double>(barrier, dE);
}
