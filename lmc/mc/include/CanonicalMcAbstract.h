/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: pravendra12                                                                 *
 * @Last Modified time: 12/01/24 10:30 AM                                                         *
 **************************************************************************************************/

/*! \file CanonicalMcAbstract.h
 *  \brief File for  CanonicalMcAbstract class declaration.
*/

#ifndef LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_
#define LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_
#include <random>
#include <mpi.h>
#include <omp.h>
#include "McAbstract.h"
#include "PotentialEnergyEstimator.h"
#include "VacancyMigrationPredictor.h"

namespace mc {

class CanonicalMcAbstract : public McAbstract {
 public:

  /*! \brief Constructor for Canonical Monte Carlo Abstract class.
   * 
   * Initializes the canonical Monte Carlo simulation with the specified 
   * configurations and parameters.
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
  CanonicalMcAbstract(Config config,
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
   * \brief Starts the simulation process.
   */
  void Simulate() override = 0;

 protected:
  // void UpdateTemperature();

  /*! 
   * \brief Dumps the current simulation state.
   */
  virtual void Dump() const;

  /*! \brief Generates a Lattice Id jump pair.
   *  \returns A pair of lattice Ids.
   */
  std::pair<size_t, size_t> GenerateLatticeIdJumpPair();
  
  /*! \brief Generates a Vacancy and Atom Jump Pair.
   *  \returns A pair of vacancy and lattice Ids.
   */
  std::pair<size_t, size_t> GenerateVacancyAtomJumpPair();

  /*! 
   * \brief Selects an event based on the lattice ID pair and energy change.
   * 
   * This function determines whether to accept or reject the proposed event, 
   * based on Metropolis Algorithm.
   * 
   * \param lattice_id_jump_pair Pair of lattice Ids.
   * \param dE                   Energy change associated with the event.
   */
  void SelectEvent(const std::pair<size_t, size_t> &lattice_id_jump_pair, 
                   double dE);
  
  /// Predictor for energy change
  const PotentialEnergyEstimator energy_change_predictor_;

  VacancyMigrationPredictor vacancy_migration_predictor_;

  /// Random Lattice Id Generator
  mutable std::uniform_int_distribution<size_t> atom_index_selector_;
};

} // mc

#endif //LMC_MC_INCLUDE_CANONICALMCABSTRACT_H_
