/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! \file CanonicalMcOmp.h
 *  @brief File for CanonicalMcOmp class declaration.
 */

#ifndef LMC_LMC_MC_INCLUDE_CANONICALMCOMP_H_
#define LMC_LMC_MC_INCLUDE_CANONICALMCOMP_H_

#include <random>
#include <mpi.h>
#include <omp.h>
#include "CanonicalMcAbstract.h"

using namespace std;

namespace mc
{
  class CanonicalMcOmp : public CanonicalMcAbstract
  {
  public:
    CanonicalMcOmp(TiledSupercell tiledSupercell,
                   unsigned long long int logDumpSteps,
                   unsigned long long int configDumpStep,
                   unsigned long long int maximumSteps,
                   unsigned long long int thermodynamicAveragingSteps,
                   unsigned long long int restartSteps,
                   double restartEnergy,
                   double temperature,
                   EnergyPredictorTLMC &energyChangePredictor);

    /** @brief Simulate canonical monte carlo simulation
     */
    void Simulate() override;

  private:
    void BuildEventVector();

    vector<pair<pair<LatticeSiteMapping, LatticeSiteMapping>, double>> event_vector_{};
    unordered_set<size_t> unavailable_position_{};

    size_t num_threads_{};
  };

} // mc

#endif // LMC_LMC_MC_INCLUDE_CANONICALMCOMP_H_
