/**************************************************************************************************
 * Copyright (c) 2023-2024. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: pravendra12                                                                 *
 * @Last Modified time: 12/01/24 6:35 PM                                                          *
 **************************************************************************************************/

/*! \file KineticMcAbstract.h
 *  \brief File for Kinetic Monte Carlo Abstract class declaration.
 */

#ifndef LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_

#include <random>
#include <omp.h>
#include <mpi.h>
#include <Eigen/Dense>
#include "McAbstract.h"
#include "JumpEvent.h"
#include "VacancyMigrationPredictor.h"
#include "TimeTemperatureInterpolator.h"
// #include "RateCorrector.hpp"


namespace mc {

class KineticMcFirstAbstract : public McAbstract {
    // it is lattice id jump pair based
  public:

   /*!
    * \brief Constructor for KineticMcFirstAbstract.
    * 
    * Initializes the kinetic Monte Carlo simulation with the given parameters.
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
    KineticMcFirstAbstract(Config config,
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

    /*!
     * \brief Destructor for KineticMcFirstAbstract.
     */
    ~KineticMcFirstAbstract() override;

    /*!
     * \brief Deleted copy constructor.
     */
    KineticMcFirstAbstract(const KineticMcFirstAbstract &) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    void operator=(const mc::KineticMcFirstAbstract &) = delete;

    /*!
     * \brief Starts the Kinetic Monte Carlo Simulation.
     */
    void Simulate() override;

  protected:
    
    /*!
     * \brief Update the temperature based on the current temperature.
     */
    void UpdateTemperature();
    // double GetTimeCorrectionFactor();

    /*!
     * \brief Dumps the current simulation state.
     */
    virtual void Dump() const;

    /*! \brief Selects an event to simulate based on rates.
     *  \return The lattice Id which will jump to vacant site.
     */
    size_t SelectEvent() const;

    /*! \brief Debugging utility to log one-step simulation time.
     *  \param one_step_time Time taken for a single simulation step.
     */
    void Debug(double one_step_time) const;

    /*! \brief Builds the event list for possible transitions.
     * 
     *  This function must be implemented by derived classes.
     */
    virtual void BuildEventList() = 0;

    /*! \brief Calculates the time increment for the simulation step.
     *         This function must be implemented by derived classes.
     * 
     *  \return Time increment.
     */
    virtual double CalculateTime() = 0;

    /*!
     * \brief Performs one step of the simulation.
     */
    virtual void OneStepSimulation();

    /// Dynamic size of the event list.
    size_t kEventListSize;

    // Helpful properties

    /// Vacancy Migration Energy Predictor.
    VacancyMigrationPredictor vacancy_migration_predictor_;

    /// Time Temperature Interpolator
    const pred::TimeTemperatureInterpolator time_temperature_interpolator_;

    /// Indicates if time-temperature interpolation is used.
    const bool is_time_temperature_interpolator_;

    /// Rate Corrector
    // const pred::RateCorrector rate_corrector_;

    /// Indicates if rate correction is used.
    const bool is_rate_corrector_;

    /// Vacancy Lattice Id.
    size_t vacancy_lattice_id_;

    /// Vacancy Trajectory Vector
    Eigen::RowVector3d vacancy_trajectory_;

    /// Jump events vector.
    std::vector<JumpEvent> event_k_i_list_{};

    /// Selected jump event.
    JumpEvent event_k_i_{};

    /// Total Rate of the events.
    double total_rate_k_{0.0}; // k would be same for all
};


/*! \brief Abstract class for Kinetic Monte Carlo Simulation.
 * 
 *  Based on Second-Order Residence Time Algorithm.
 */
class KineticMcChainAbstract : public KineticMcFirstAbstract {
  public:
    
    /*!
    * \brief Constructor for KineticMcChainAbstract.
    * 
    * Initializes the abstract kinetic Monte Carlo simulation based on 
    * Second-Order Residence Time Algorithm with the given parameters.
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
    KineticMcChainAbstract(Config config,
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
    
    /*!
     * \brief Destructor for KineticMcChainAbstract.
     */
    ~KineticMcChainAbstract() override;

    /*!
     * \brief Deleted copy constructor.
     */
    KineticMcChainAbstract(const KineticMcChainAbstract &) = delete;

    /*!
     * \brief Deleted assignment operator.
     */
    void operator=(const mc::KineticMcChainAbstract &) = delete;

  protected:
    
    /*!
     * @brief Performs one step of the simulation.
     */
    void OneStepSimulation() override;

    // helpful properties

    /// Lattice ID of the previous jump site.
    size_t previous_j_lattice_id_;

    /// Total rate for jump to nearest neighbours.
    double total_rate_i_{0.0}; // i would be different

    /// Neighbours of choosen lattice Id for the jump.
    ///  j -> k -> i ->l
    ///  k : current position
    std::vector<size_t> l_lattice_id_list_{};

    MPI_Op mpi_op_{};
    MPI_Datatype mpi_datatype_{};
};

/*!
 * \struct MpiData
 * \brief Structure for storing data used in MPI custom reduction.
 */
struct MpiData {
  double beta_bar_k{0.0};
  double beta_k{0.0};
  double gamma_bar_k_j{0.0};
  double gamma_k_j{0.0};
  double beta_k_j{0.0};
  double alpha_k_j{0.0};
  double ts_numerator{0.0};
  double ts_j_numerator{0.0};
};

/*!
 * \brief Performs a custom reduction operation for MpiData.
 *
 * \param input_buffer Input buffer.
 * \param output_buffer Output buffer.
 * \param len Number of elements.
 * \param datatype MPI datatype.
 */
inline void DataSum(void *input_buffer,
                    void *output_buffer,
                    int *len,
                    [[maybe_unused]] MPI_Datatype *datatype) {
  auto *input = static_cast<MpiData *>(input_buffer);
  auto *output = static_cast<MpiData *>(output_buffer);
  for (int i = 0; i < *len; ++i) {
    output[i].beta_bar_k += input[i].beta_bar_k;
    output[i].beta_k += input[i].beta_k;
    output[i].gamma_bar_k_j += input[i].gamma_bar_k_j;
    output[i].gamma_k_j += input[i].gamma_k_j;
    output[i].beta_k_j += input[i].beta_k_j;
    output[i].alpha_k_j += input[i].alpha_k_j;
    output[i].ts_numerator += input[i].ts_numerator;
    output[i].ts_j_numerator += input[i].ts_j_numerator;
  }
}


inline void DefineStruct(MPI_Datatype *datatype) {
  const int count = 8;
  int block_lens[count];
  MPI_Datatype types[count];
  MPI_Aint displacements[count];

  for (int i = 0; i < count; i++) {
    types[i] = MPI_DOUBLE;
    block_lens[i] = 1;
  }
  displacements[0] = offsetof(MpiData, beta_bar_k);
  displacements[1] = offsetof(MpiData, beta_k);
  displacements[2] = offsetof(MpiData, gamma_bar_k_j);
  displacements[3] = offsetof(MpiData, gamma_k_j);
  displacements[4] = offsetof(MpiData, beta_k_j);
  displacements[5] = offsetof(MpiData, alpha_k_j);
  displacements[6] = offsetof(MpiData, ts_numerator);
  displacements[7] = offsetof(MpiData, ts_j_numerator);
  MPI_Type_create_struct(count, block_lens, displacements, types, datatype);
  MPI_Type_commit(datatype);
}
} // mc

#endif //LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_
