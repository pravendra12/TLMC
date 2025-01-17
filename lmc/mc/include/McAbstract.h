/**************************************************************************************************
 * Copyright (c) 2023-2024. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: pravendra12                                                                 *
 * @Last Modified time: 12/02/24 9:26 AM                                                          *
 **************************************************************************************************/

/*! \file McAbstract.h
 *  \brief File for Monte Carlo Abstract class declaration.
 */

#ifndef LMC_MC_INCLUDE_MCABSTRACT_H_
#define LMC_MC_INCLUDE_MCABSTRACT_H_
#include <random>
#include <omp.h>
#include <mpi.h>
#include "Config.h"
#include "ThermodynamicAveraging.h"

class McAbstract {
 public:
  
  /*! \brief Constructor for Monte Carlo Abstract class.
   *
   *  \param config                        Configuration used for simulation.
   *  \param supercell_config              Configuration used for training the 
   *                                       Cluster Expansion Model.
   *  \param log_dump_steps                Number of steps between logging the 
   *                                       simulation progress.
   *  \param config_dump_steps             Number of steps between dumping the 
   *                                       configuration.
   *  \param maximum_steps                 Maximum number of steps to perform in 
   *                                       the simulation.
   *  \param thermodynamic_averaging_steps Number of steps for thermodynamic 
   *                                       averaging.
   *  \param restart_steps                 Step count for restarting the 
   *                                       simulation.
   *  \param restart_time                  Time value for restarting the 
   *                                       simulation.
   *  \param restart_energy                Energy value for restarting the 
   *                                       simulation.
   *  \param temperature                   Constant temperature for the 
   *                                       simulation (in Kelvin).
   *  \param element_set                   Set of elements used in the 
   *                                       configuration.
   *  \param max_cluster_size              Maximum size of clusters to be 
   *                                       considered in the simulation.
   *  \param max_bond_order                Maximum bond order used for 
   *                                       determining clusters.
   *  \param json_coefficients_filename    Path to the JSON file containing 
   *                                       cluster interaction coefficients.
   *  \param log_filename                  Name of log file.
   */
  McAbstract(Config config,
             Config supercell_config,
             unsigned long long int log_dump_steps,
             unsigned long long int config_dump_steps,
             unsigned long long int maximum_steps,
             unsigned long long int thermodynamic_averaging_steps,
             unsigned long long int restart_steps,
             double restart_energy,
             double restart_time,
             double temperature,
             const std::set<Element> &element_set,
             const size_t max_cluster_size,
             const size_t max_bond_order,
             const std::string &json_coefficients_filename,
             const std::string &log_filename);

  /*!
   *  \brief Destructor for McAbstract.
   */
  virtual ~McAbstract();

  /*!
   *  \brief Deleted copy constructor to prevent copying.
   */
  McAbstract(const McAbstract &) = delete;

  /**
   *  \brief Deleted assignment operator to prevent copying.
   */
  void operator=(const McAbstract &) = delete;

  /*!
   *  \brief Starts the simulation.
   */
  virtual void Simulate() = 0;

 protected:

  /// Store the configuration.
  Config config_;

  // simulation parameters

  /// Number of steps between logging the simulation process.
  const unsigned long long int log_dump_steps_;

  /// Number of steps between dumping the configuration.
  const unsigned long long int config_dump_steps_;

  /// Maximum number of steps for simulation.
  const unsigned long long int maximum_steps_;

  // simulation statistics

  unsigned long long int steps_;
  double energy_;
  double absolute_energy_;
  double time_;
  double temperature_;
  double beta_;
  mutable bool is_restarted_;

  // helpful properties
  // mc::ThermodynamicAveraging thermodynamic_averaging_;
  mutable std::mt19937_64 generator_;
  mutable std::uniform_real_distribution<double> unit_distribution_;
  mutable std::ofstream ofs_;

  int world_rank_{-1};
  int world_size_{-1};
};

#endif //LMC_MC_INCLUDE_MCABSTRACT_H_
