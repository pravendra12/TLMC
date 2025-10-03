#ifndef LMC_PRED_INCLUDE_ENERGYPREDICTORTLMC_H_
#define LMC_PRED_INCLUDE_ENERGYPREDICTORTLMC_H_

#include <chrono>
#include "Config.h"
#include "PrintUtility.h"
#include "ClusterExpansionParameters.h"
#include "SymmetricCEPredictorTLMC.h"
#include "LVFEPredictorTLMC.h"
#include "TiledSupercell.h"

using namespace std;

class EnergyPredictorTLMC
{
public:
  EnergyPredictorTLMC(
      SymmetricCEPredictorTLMC &symCEEnergyPredictor,
      LVFEPredictorTLMC &lvfePredictor);

  // Applied for any pair
  double GetEnergyChange(
    const TiledSupercell &tiledSupercell, 
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair);

  // Specifically for <vacancy, atom> pair);
  double GetEnergyChangeWithVacancy(
      const TiledSupercell &tiledSupercell,
      const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair);

  void ProfileEnergyChange(
    const TiledSupercell &tiledSupercell, 
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair, 
    int numThreads);

private:
  SymmetricCEPredictorTLMC &symCEEnergyPredictor_;
  LVFEPredictorTLMC &lvfePredictor_;

  const vector<string> &allowedElements_;
};

#endif // LMC_PRED_INCLUDE_ENERGYPREDICTORTLMC_H_