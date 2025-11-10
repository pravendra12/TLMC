/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/**
 * @file KineticMcChainOmpi.h
 * @brief Declaration of the KineticMcChainOmpi class for OMPI-based kinetic Monte Carlo simulations.
 *        j -> k -> i -> l
 *             |
 *        current position
 */

#ifndef LMC_LMC_MC_INCLUDE_KINETICMCCHAINOMPI_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCCHAINOMPI_H_

#include <random>
#include <omp.h>
#include <mpi.h>
#include "JumpEvent.h"
#include "KineticMcAbstract.h"

using namespace std;

namespace mc
{

  class KineticMcChainOmpi : public KineticMcChainAbstract
  {
  public:
    /**
     * @brief Constructor for KineticMcFirstOmpi.
     *
     * Sets up the kinetic Monte Carlo simulation using the Second-Order Residence Time Algorithm.
     *
     */
    KineticMcChainOmpi(TiledSupercell tiledSupercell,
                       const unsigned long long int logDumpSteps,
                       const unsigned long long int configDumpSteps,
                       const unsigned long long int maximumSteps,
                       const unsigned long long int thermodynamicAveragingSteps,
                       const unsigned long long int restartSteps,
                       const double restartEnergy,
                       const double restartTime,
                       const double temperature,
                       VacancyMigrationPredictorTLMC &vacancyMigrationPredictor,
                       const string &timeTemperatureFilename,
                       unique_ptr<RateCorrector> &rateCorrector,
                       const Eigen::RowVector3d &vacancyTrajectory);

  protected:
    /*!
     * \brief Builds the event list for the simulation.
     *
     * Implements the mechanism to populate the list of possible transitions
     * and their associated rates based on the current simulation state.
     */
    void BuildEventList() override;

    /*!
     * \brief Calculates the time increment for the simulation step.
     *
     * Computes the time increment based on the total transition rate
     * and updates the simulation time accordingly.
     *
     * \return Time increment for the current simulation step.
     */
    double CalculateTime() override;
  };

} // mc

#endif // LMC_LMC_MC_INCLUDE_KINETICMCCHAINOMPI_H_
