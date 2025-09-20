/*******************************************************************************
 * Copyright (c) 2025. All rights reserved.
 * @Author: Pravendra Patel
 * @Date:    2025-06-02
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-02
 *******************************************************************************/

/*! @file VacancyMigrationPredictorTLMC.h
    @brief File contains declaration of vacancy migration predictor class.
*/

#ifndef LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORTLMC_H_
#define LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORTLMC_H_

#include "KRAPredictorTLMC.h"
#include "EnergyPredictorTLMC.h"
#include "TiledSupercell.h"
#include "LatticeSiteMapping.hpp"

using namespace std;

class VacancyMigrationPredictorTLMC
{
public:

  VacancyMigrationPredictorTLMC(
      KRAPredictorTLMC &eKRAPredictor, 
      EnergyPredictorTLMC &energyPredictor);
  
  /*! @brief Returns barrier and energy change due to vacancy migration
      @param config Configuration for which barrier and energy change to computed
      @param latticeIdJumpPair Lattice Id Jump Pair
  */
  pair<double, double> GetBarrierAndDeltaE(
      const TiledSupercell &tiledSupercell,
      const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeIdJumpPair);

private:
  /*! @brief KRA predictor object
   */
  KRAPredictorTLMC &eKRAPredictor_;

  /*! @brief Energy change predictor object
   */
  EnergyPredictorTLMC &energyPredictor_;
};

#endif // LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORTLMC_H_
