/*******************************************************************************
 * Copyright (c) 2025. All rights reserved.
 * @Author: Pravendra Patel
 * @Date:    2025-06-02
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-02
 *******************************************************************************/

/*! @file VacancyMigrationPredictor.h
    @brief File contains implementation of vacancy migration predictor class.
*/

#include "VacancyMigrationPredictor.h"

VacancyMigrationPredictor::VacancyMigrationPredictor(
    const ClusterExpansionParameters &ceParams,
    const Config &config) : eKRAPredictor_(ceParams,
                                           config),
                            energyPredictor_(
                                ceParams,
                                config)
{
}

pair<double, double> VacancyMigrationPredictor::GetBarrierAndDeltaE(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair) 
{
  double eKRA = eKRAPredictor_.GetKRA(config, latticeIdJumpPair);

  double dE = energyPredictor_.GetDeMigration(
      config,
      latticeIdJumpPair);

  double barrier = eKRA + (dE / 2);

  return pair<double, double>(barrier, dE);
}
