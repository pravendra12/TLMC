/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/**
 * @file KineticMcFirstMpi.h
 * @brief Declaration of the KineticMcFirstMpi class for MPI-based kinetic Monte Carlo simulations.
 */

#ifndef LMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_
#define LMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_

#include "McAbstract.h"
#include "KineticMcAbstract.h"
#include "VacancyMigrationPredictorTLMC.h"

using namespace std;

namespace mc
{

  class KineticMcFirstMpi : public KineticMcFirstAbstract
  {
  public:
    /**
     * @brief Constructor for KineticMcFirstMpi.
     *
     * Sets up the kinetic Monte Carlo simulation using the First-Order Residence Time Algorithm.
     */
    KineticMcFirstMpi(TiledSupercell tiledSupercell,
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

    /*!
     * @brief Deconstructor for KineticMcFirstMpi class.
     */
    ~KineticMcFirstMpi() override;

  protected:
    /*!
     * @brief Builds the event list for the simulation.
     *
     * Implements the mechanism to populate the list of possible transitions
     * and their associated rates based on the current simulation state.
     */
    void BuildEventList() override;

    /*!
     * @brief Calculates the time increment for the simulation step.
     *
     * Computes the time increment based on the total transition rate
     * and updates the simulation time accordingly.
     *
     * @return Time increment for the current simulation step.
     */
    double CalculateTime() override;
  };
} // mc

#endif // LMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_