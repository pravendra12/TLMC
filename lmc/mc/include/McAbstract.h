/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/**
 * @file McAbstract.h
 * @brief File contains declaration of McAbstract Class.
 */

#ifndef LMC_MC_INCLUDE_MCABSTRACT_H_
#define LMC_MC_INCLUDE_MCABSTRACT_H_

#include <random>
#include <utility>
#include <omp.h>
#include <mpi.h>
#include "TiledSupercell.h"
#include "ThermodynamicAveraging.h"

using namespace std;

class McAbstract
{
public:

  McAbstract(TiledSupercell tiledSupercell,
             unsigned long long int logDumpSteps,
             unsigned long long int configDumpSteps,
             unsigned long long int maximumSteps,
             unsigned long long int thermodynamicAveragingSteps,
             unsigned long long int restartSteps,
             double restartEnergy,
             double restartTime,
             double temperature,
             const string &logFilename);

  /**
   *  \brief Destructor for McAbstract.
   */
  virtual ~McAbstract();

  /**
   *  \brief Deleted copy constructor to prevent copying.
   */
  McAbstract(const McAbstract &) = delete;

  /**
   *  \brief Deleted assignment operator to prevent copying.
   */
  void operator=(const McAbstract &) = delete;

  /**
   *  \brief Starts the simulation.
   */
  virtual void Simulate() = 0;

protected:
  /**
   * @brief Stores the configuration.
   *
   */
  TiledSupercell tiledSupercell_;

  // simulation parameters

  /**
   * @brief Number of steps between logging the simulation process.
   *
   */
  const unsigned long long int logDumpSteps_;

  /**
   * @brief Number of steps between dumping the configuration.
   *
   */
  const unsigned long long int configDumpSteps_;

  /**
   * @brief Maximum number of steps for simulation.
   *
   */
  const unsigned long long int maximumSteps_;

  /**
   * @brief simulation statistics
   *
   */

  unsigned long long int steps_;
  double energy_;
  double absolute_energy_;
  double time_;
  double temperature_;
  double beta_;
  mutable bool is_restarted_;

  mc::ThermodynamicAveraging thermodynamicAveraging_;
  mutable mt19937_64 generator_;
  mutable uniform_real_distribution<double> unitDistribution_;
  mutable ofstream ofs_;

  int world_rank_{-1};
  int world_size_{-1};
};

#endif // LMC_MC_INCLUDE_MCABSTRACT_H_
