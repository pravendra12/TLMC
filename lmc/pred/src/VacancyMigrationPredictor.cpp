#include "VacancyMigrationPredictor.h"

VacancyMigrationPredictor::VacancyMigrationPredictor(
    const string &predictorFilename,
    const Config &referenceConfig,
    const Config &supercellConfig,
    const set<Element> &elementSet,
    size_t maxClusterSize,
    size_t maxBondOrder) : barrierPredictor_(referenceConfig,
                                             elementSet,
                                             predictorFilename),
                           energyChangePredictor_(predictorFilename,
                                                  referenceConfig,
                                                  supercellConfig,
                                                  elementSet,
                                                  maxClusterSize,
                                                  maxBondOrder)
{}

pair<double, double> VacancyMigrationPredictor::getBarrierAndEnergyChange(
    const Config &config, 
    const pair<size_t, size_t> &latticeIdJumpPair) const
{
  
  double barrier = barrierPredictor_.GetBarrier(config, 
                                                latticeIdJumpPair);
                                               
  double dE = energyChangePredictor_.GetDeThreadSafe(config,
                                                     latticeIdJumpPair);
 
  return {barrier, dE};
}
