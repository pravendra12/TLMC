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
#include "PotentialEnergyEstimator.h"

namespace mc
{
  class CanonicalMcSerial : public CanonicalMcAbstract
  {
  public:
    /**
     * @brief Constructor for the Canonical Monte Carlo Omp class.
     *
     * @param config Configuration used for the simulation.
     * @param supercellConfig Configuration used to train the Cluster Expansion model.
     * @param logDumpSteps Number of steps between logging simulation progress.
     * @param configDumpSteps Number of steps between dumping the configuration.
     * @param maximumSteps Maximum number of Monte Carlo steps to run.
     * @param thermodynamicAveragingSteps Number of steps to use for thermodynamic averaging.
     * @param restartSteps Step index for restarting the simulation.
     * @param restartEnergy Energy value to resume the simulation from.
     * @param temperature Constant simulation temperature (in Kelvin).
     * @param elementSet Set of elements present in the configuration.
     * @param predictorFilename Path to the JSON file containing cluster interaction coefficients.
     */

    CanonicalMcSerial(Config config,
                      Config supercellConfig,
                      unsigned long long int logDumpSteps,
                      unsigned long long int configDumpStep,
                      unsigned long long int maximumSteps,
                      unsigned long long int thermodynamicAveragingSteps,
                      unsigned long long int restartSteps,
                      double restartEnergy,
                      double temperature,
                    const ClusterExpansionParameters &ceParams);

    /** @brief Simulate canonical monte carlo simulation
     */
    void Simulate() override;
  };
} // mc

#endif // LMC_MC_INCLUDE_CANONICALMCSERIAL_H_
