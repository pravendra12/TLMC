/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! \file CanonicalMcAbstract.h
 *  @brief File for CanonicalMcAbstract class declaration.
 */

#ifndef LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_
#define LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_

#include <random>
#include <mpi.h>
#include <omp.h>
#include <utility>
#include <chrono>
#include "McAbstract.h"
#include "EnergyPredictorTLMC.h"

using namespace std;

namespace mc
{
  class CanonicalMcAbstract : public McAbstract
  {
  public:
    /**
     * @brief Constructor for the Canonical Monte Carlo abstract class.
     *
     * Initializes the Canonical Monte Carlo simulation with the given configuration,
     * training setup, and simulation parameters.
     *
     */

    CanonicalMcAbstract(TiledSupercell tiledSupercell,
                        unsigned long long int logDumpSteps,
                        unsigned long long int configDumpSteps,
                        unsigned long long int maximumSteps,
                        unsigned long long int thermodynamicAveragingSteps,
                        unsigned long long int restartSteps,
                        double restartEnergy,
                        double temperature,
                        EnergyPredictorTLMC &energyChangePredictor);

    /*!
     * @brief Starts the simulation process.
     */
    void Simulate() override = 0;

  protected:
    // void UpdateTemperature();

    /*!
     * @brief Dumps the current simulation state.
     */
    virtual void Dump() const;

    /*! @brief Generates a lattice Id jump pair.
     *  \returns A pair of lattice Ids.
     */
    pair<LatticeSiteMapping, LatticeSiteMapping> GenerateLatticeSiteIdJumpPair();

    pair<LatticeSiteMapping, LatticeSiteMapping> GenerateVacancyLatticeSiteIdJumpPair();

    /*!
     * @brief Selects an event based on the lattice ID pair and energy change.
     *
     * This function determines whether to accept or reject the proposed event,
     * based on Metropolis Algorithm.
     *
     * \param lattice_id_jump_pair Pair of lattice Ids.
     * \param dE                   Energy change associated with the event.
     */
    void SelectEvent(const pair<LatticeSiteMapping, LatticeSiteMapping> &lattice_id_jump_pair,
                     const double dE);


    EnergyPredictorTLMC &energyChangePredictor_;

    /** @brief Random Lattice Id Generator
     */
    mutable uniform_int_distribution<size_t> atomIndexSelector_;
  };
}

#endif // LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_
