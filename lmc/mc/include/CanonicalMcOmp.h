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
#include "PotentialEnergyEstimator.h"

using namespace std;

namespace mc
{
  class CanonicalMcOmp : public CanonicalMcAbstract
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
    CanonicalMcOmp(Config config,
                   Config supercellConfig,
                   unsigned long long int logDumpSteps,
                   unsigned long long int configDumpStep,
                   unsigned long long int maximumSteps,
                   unsigned long long int thermodynamicAveragingSteps,
                   unsigned long long int restartSteps,
                   double restartEnergy,
                   double temperature,
                   const set<Element> &elementSet,
                   const string &predictorFilename);

    /** @brief Simulate canonical monte carlo simulation
     */
    void Simulate() override;

  private:
    void BuildEventVector();

    vector<pair<pair<size_t, size_t>, double>> event_vector_{};
    unordered_set<size_t> unavailable_position_{};

    size_t num_threads_{};
  };

} // mc

#endif // LMC_LMC_MC_INCLUDE_CANONICALMCOMP_H_
