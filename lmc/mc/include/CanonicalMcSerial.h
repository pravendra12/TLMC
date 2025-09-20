/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! \file CanonicalMcSerial.h
 *  @brief File for CanonicalMcSerial class implementation.
 */

#ifndef LMC_MC_INCLUDE_CANONICALMCSERIAL_H_
#define LMC_MC_INCLUDE_CANONICALMCSERIAL_H_

#include <random>
#include <mpi.h>
#include <omp.h>
#include "CanonicalMcAbstract.h"

namespace mc
{
  class CanonicalMcSerial : public CanonicalMcAbstract
  {
  public:
    /**
     * @brief Constructor for the Canonical Monte Carlo Omp class.
     *
     */

    CanonicalMcSerial(TiledSupercell tiledSupercell,
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
  };
} // mc

#endif // LMC_MC_INCLUDE_CANONICALMCSERIAL_H_
