/*******************************************************************************
 * Copyright (c) 2025. All rights reserved.
 * @Author: Pravendra Patel
 * @Date:    2025-06-02
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-02
 *******************************************************************************/

/*! @file VacancyMigrationPredictorTLMC.h
    @brief File contains implementation of vacancy migration predictor class.
*/

#include "VacancyMigrationPredictorTLMC.h"

VacancyMigrationPredictorTLMC::VacancyMigrationPredictorTLMC(
    KRAPredictorTLMC &eKRAPredictor,
    EnergyPredictorTLMC &energyPredictor) : eKRAPredictor_(eKRAPredictor),
                                            energyPredictor_(energyPredictor)
{
}

pair<double, double> VacancyMigrationPredictorTLMC::GetBarrierAndDeltaE(
    const TiledSupercell &tiledSupercell,
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair)
{
  double eKRA = eKRAPredictor_.GetKRA(tiledSupercell, latticeSiteJumpPair);

  double dE = energyPredictor_.GetEnergyChange(
      tiledSupercell,
      latticeSiteJumpPair);

  double barrier = eKRA + (dE / 2);

  return pair<double, double>(barrier, dE);
}
