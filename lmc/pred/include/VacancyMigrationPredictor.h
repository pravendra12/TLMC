/*******************************************************************************
 * Copyright (c) 2025. All rights reserved.
 * @Author: Pravendra Patel
 * @Date:    2025-06-02
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-02
 *******************************************************************************/

/*! @file VacancyMigrationPredictor.h
    @brief File contains declaration of vacancy migration predictor class.
*/

#ifndef LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_
#define LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_

#include "KRAPredictor.h"
#include "PotentialEnergyEstimator.h"
#include "ClusterExpansionParameters.h"

using namespace std;

class VacancyMigrationPredictor
{
public:

  VacancyMigrationPredictor(
      const ClusterExpansionParameters &ceParams, 
      const Config &config);
  
  /*! @brief Returns barrier and energy change due to vacancy migration
      @param config Configuration for which barrier and energy change to computed
      @param latticeIdJumpPair Lattice Id Jump Pair
  */
  pair<double, double> GetBarrierAndDeltaE(
      const Config &config,
      const pair<size_t, size_t> &latticeIdJumpPair);

private:
  /*! @brief KRA predictor object
   */
  KRAPredictor eKRAPredictor_;

  /*! @brief Energy change predictor object
   */
  PotentialEnergyEstimator energyPredictor_;
};

#endif // LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_
