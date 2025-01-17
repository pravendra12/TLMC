/**************************************************************************************************
 * Copyright (c) 2023-2024. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: pravendra12                                                                 *
 * @Last Modified time: 12/02/24 9:26 AM                                                          *
 **************************************************************************************************/

/*! \file KineticMcChainOmp.h
 *  \brief File for Kinetic Monte Carlo class declaration.
 */

#ifndef LMC_LMC_MC_INCLUDE_KINETICMCCHAINOMPI_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCCHAINOMPI_H_
#include <random>
#include <omp.h>
#include <mpi.h>
#include "VacancyMigrationPredictor.h"
#include "JumpEvent.h"
#include "KineticMcAbstract.h"
namespace mc {

//  j -> k -> i -> l
//       |
// current position

class KineticMcChainOmpi : public KineticMcChainAbstract {
  public:

    /*!
    * \brief Constructor for KineticMcChainOmpi
    * 
    * Initializes the kinetic Monte Carlo simulation based on Second-Order 
    * Residence Time Algorithm with the given parameters.
    * 
    * \param config                        Simulation configuration.
    * \param supercell_config              Training configuration for the 
    *                                      Cluster Expansion Model.
    * \param log_dump_steps                Steps between logging progress.
    * \param config_dump_steps             Steps between configuration dumps.
    * \param maximum_steps                 Maximum simulation steps.
    * \param thermodynamic_averaging_steps Steps for thermodynamic averaging.
    * \param restart_steps                 Steps for restarting the simulation.
    * \param restart_energy                Restart energy.
    * \param restart_time                  Restart time.
    * \param temperature                   Simulation temperature (in Kelvin).
    * \param element_set                   Set of elements involved in the 
    *                                      simulation.
    * \param max_cluster_size              Maximum size of clusters to consider.
    * \param max_bond_order                Maximum bond order for determining 
    *                                      clusters.
    * \param json_coefficients_filename    Path to JSON file with cluster 
    *                                      interaction coefficients.
    * \param time_temperature_filename     Path to time-temperature data file.
    * \param is_rate_corrector             Whether rate correction need to
    *                                      applied.
    * \param vacancy_trajectory            Initial vacancy trajectory vector.
    */
    KineticMcChainOmpi(Config config,
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
                       const std::string &time_temperature_filename,
                       bool is_rate_corrector,
                       const Eigen::RowVector3d &vacancy_trajectory);

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

#endif //LMC_LMC_MC_INCLUDE_KINETICMCCHAINOMPI_H_
