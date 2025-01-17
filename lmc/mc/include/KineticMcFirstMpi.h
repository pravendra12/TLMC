/**************************************************************************************************
 * Copyright (c) 2020-2024. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 11/26/24 4:35 PM                                                                        *
 * @Last Modified by: pravendra                                                                   *
 * @Last Modified time: 11/26/24 4:35 PM                                                          *
 **************************************************************************************************/

/*! \file  KineticMcFirstMpi.h
 *  \brief File for the KineticMcFirstMpi class declaration.
 */

#ifndef LMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_
#define LMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_

#include <random>
#include <omp.h>
#include <mpi.h>
#include "JumpEvent.h"
#include "KineticMcAbstract.h"
#include "VacancyMigrationPredictor.h"

namespace mc {

class KineticMcFirstMpi : public KineticMcFirstAbstract {
  public:

    /*!
    * \brief Constructor for KineticMcFirstMpi
    * 
    * Initializes the kinetic Monte Carlo simulation based on First-Order 
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
    KineticMcFirstMpi(Config config,
                     Config supercell_config,
                     const unsigned long long int log_dump_steps,
                     const unsigned long long int config_dump_steps,
                     const unsigned long long int maximum_steps,
                     const unsigned long long int thermodynamic_averaging_steps,
                     const unsigned long long int restart_steps,
                     const double restart_energy,
                     const double restart_time,
                     const double temperature,
                     const std::set<Element> &element_set,
                     const size_t max_cluster_size,
                     const size_t max_bond_order,
                     const std::string &json_coefficients_filename,
                     const std::string &time_temperature_filename,
                     const bool is_rate_corrector,
                     const Eigen::RowVector3d &vacancy_trajectory);

    /*!
     * \brief Deconstructor for KineticMcFirstMpi class.
     */
    ~KineticMcFirstMpi() override;

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

#endif //LMC_MC_INCLUDE_KINETICMCFIRSTMPI_H_