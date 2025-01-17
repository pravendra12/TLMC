/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: pravendra12                                                                 *
 * @Last Modified time: 12/01/24 10:35 AM                                                         *
 **************************************************************************************************/


/*! \file CanonicalMcSerial.h
 *  \brief File for Canonical Monte Carlo Serial class declaration.
 */

#ifndef LMC_MC_INCLUDE_CANONICALMCSERIAL_H_
#define LMC_MC_INCLUDE_CANONICALMCSERIAL_H_
#include <random>
#include <mpi.h>
#include <omp.h>
#include "CanonicalMcAbstract.h"
#include "PotentialEnergyEstimator.h"

namespace mc {
class CanonicalMcSerial : public CanonicalMcAbstract {
 public:

  /*! \brief Constructor for Canonical Monte Carlo Serial Class
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
   */
  CanonicalMcSerial(Config config,
                    Config supercell_config,
                    unsigned long long int log_dump_steps,
                    unsigned long long int config_dump_steps,
                    unsigned long long int maximum_steps,
                    unsigned long long int thermodynamic_averaging_steps,
                    unsigned long long int restart_steps,
                    double restart_energy,
                    double temperature,
                    const std::set<Element> &element_set,                                 
                    const size_t max_cluster_size,
                    const size_t max_bond_order,
                    const std::string &json_coefficients_filename);

  /*!
   * \brief Starts the simulation.
   */
  void Simulate() override;
};
} // mc

#endif //LMC_MC_INCLUDE_CANONICALMCSERIAL_H_
