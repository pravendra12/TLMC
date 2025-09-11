/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/**
 * @file KineticMcAbstract.h
 * @brief File contains declaration of KineticMc abstract class.
 *
 */

#ifndef LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_

#include <random>
#include <omp.h>
#include <mpi.h>
#include <Eigen/Dense>
#include "McAbstract.h"
#include "JumpEvent.h"
#include "TimeTemperatureInterpolator.h"
#include "VacancyMigrationPredictor.h"

using namespace std;

namespace mc
{

  /**
   * @brief Abstract class for Kinetic Monte Carlo Simulation.
   *
   *  Based on First-Order Residence Time Algorithm.
   */
  class KineticMcFirstAbstract : public McAbstract
  {
  public:
    /**
     * Constructor for KineticMcFirstAbstract.
     *
     * Initializes the kinetic Monte Carlo simulation with the given parameters.
     *
     * @param config Simulation configuration.
     * @param supercellConfig Training configuration used for the Cluster Expansion Model.
     * @param logDumpSteps Steps between logging progress.
     * @param configDumpSteps Steps between configuration dumps.
     * @param maximumSteps Maximum simulation steps.
     * @param thermodynamicAveragingSteps Steps for thermodynamic averaging.
     * @param restartSteps Steps for restarting the simulation.
     * @param restartEnergy Restart energy.
     * @param restartTime Restart time.
     * @param temperature Simulation temperature (in Kelvin).
     * @param elementSet Set of elements involved in the simulation.
     * @param predictorFilename Path to JSON file with cluster interaction coefficients.
     * @param timeTemperatureFilename Path to time-temperature data file.
     * @param isRateCorrector Whether rate correction needs to be applied.
     * @param vacancyTrajectory Initial vacancy trajectory vector.
     */
    KineticMcFirstAbstract(Config config,
                           unsigned long long int logDumpSteps,
                           unsigned long long int configDumpSteps,
                           unsigned long long int maximumSteps,
                           unsigned long long int thermodynamicAveragingSteps,
                           unsigned long long int restartSteps,
                           double restartEnergy,
                           double restartTime,
                           double temperature,
                           VacancyMigrationPredictor &vacancyMigrationPredictor,
                           const string &timeTemperatureFilename,
                           bool isRateCorrector,
                           const Eigen::RowVector3d &vacancyTrajectory);

    /**
     * @brief Destructor for KineticMcFirstAbstract.
     */
    ~KineticMcFirstAbstract() override;

    /**
     * @brief Deleted copy constructor.
     */
    KineticMcFirstAbstract(const KineticMcFirstAbstract &) = delete;

    /**
     * @brief Deleted assignment operator.
     */
    void operator=(const mc::KineticMcFirstAbstract &) = delete;

    /**
     * @brief Starts the Kinetic Monte Carlo Simulation.
     */
    void Simulate() override;

  protected:
    /**
     * @brief Update the temperature based on the current temperature.
     */
    void UpdateTemperature();

    /**
     * @brief Dumps the current simulation state.
     */
    virtual void Dump() const;

    /**
     * @brief Selects an event to simulate based on rates.
     * @return The lattice Id which will jump to vacant site.
     */
    size_t SelectEvent() const;

    /** @brief Debugging utility to log one-step simulation time.
     *  @param one_step_time Time taken for a single simulation step.
     */
    void Debug(double one_step_time) const;

    /** @brief Builds the event list for possible transitions.
     *
     *  This function must be implemented by derived classes.
     */
    virtual void BuildEventList() = 0;

    /** @brief Calculates the time increment for the simulation step.
     *         This function must be implemented by derived classes.
     *
     *  @return Time increment.
     */
    virtual double CalculateTime() = 0;

    /**
     * @brief Performs one step of the simulation.
     */
    virtual void OneStepSimulation();

    /**
     * @brief  Dynamic size of the event list.
     */
    size_t kEventListSize_;

    // Helpful properties

    /**
     * @brief Vacancy Migration Energy Predictor
     */
    VacancyMigrationPredictor &vacancyMigrationPredictor_;

    /**
     * @brief Time Temperature Interpolator
     *
     */
    const TimeTemperatureInterpolator timeTemperatureInterpolator_;

    /**
     * @brief Indicates if time-temperature interpolation is used
     *
     */
    const bool isTimeTemperatureInterpolator_;

    // @brief  Rate Corrector
    // const pred::RateCorrector rate_corrector_;

    /**
     * @brief Indicates if rate correction is used.
     *
     */
    const bool isRateCorrector_;

    /**
     * @brief Vacancy Lattice Id.
     *
     */
    size_t vacancyLatticeId_;

    /**
     * @brief Vacancy Trajectory Vector
     *
     */
    Eigen::RowVector3d vacancyTrajectory_;

    /**
     * @brief Jump events vector.
     *
     */
    vector<JumpEvent> event_k_i_list_{};

    /**
     * @brief Selected jump event.
     *
     */
    JumpEvent event_k_i_{};

    /**
     * @brief Total Rate of the events.
     *
     */
    double total_rate_k_{0.0};
  };

  /**
   * @brief Abstract class for Kinetic Monte Carlo Simulation.
   *
   *  Based on Second-Order Residence Time Algorithm.
   */
  class KineticMcChainAbstract : public KineticMcFirstAbstract
  {
  public:
    /**
     * Constructor for KineticMcFirstAbstract.
     *
     * Initializes the kinetic Monte Carlo simulation with the given parameters.
     *
     * @param config Simulation configuration.
     * @param supercellConfig Training configuration used for the Cluster Expansion Model.
     * @param logDumpSteps Steps between logging progress.
     * @param configDumpSteps Steps between configuration dumps.
     * @param maximumSteps Maximum simulation steps.
     * @param thermodynamicAveragingSteps Steps for thermodynamic averaging.
     * @param restartSteps Steps for restarting the simulation.
     * @param restartEnergy Restart energy.
     * @param restartTime Restart time.
     * @param temperature Simulation temperature (in Kelvin).
     * @param elementSet Set of elements involved in the simulation.
     * @param predictorFilename Path to JSON file with cluster interaction coefficients.
     * @param timeTemperatureFilename Path to time-temperature data file.
     * @param isRateCorrector Whether rate correction needs to be applied.
     * @param vacancyTrajectory Initial vacancy trajectory vector.
     */
    KineticMcChainAbstract(Config config,
                           unsigned long long int logDumpSteps,
                           unsigned long long int configDumpSteps,
                           unsigned long long int maximumSteps,
                           unsigned long long int thermodynamicAveragingSteps,
                           unsigned long long int restartSteps,
                           double restartEnergy,
                           double restartTime,
                           double temperature,
                           VacancyMigrationPredictor &vacancyMigrationPredictor,
                           const string &timeTemperatureFilename,
                           bool isRateCorrector,
                           const Eigen::RowVector3d &vacancyTrajectory);

    /**
     * @brief Destructor for KineticMcChainAbstract.
     */
    ~KineticMcChainAbstract() override;

    /**
     * @brief Deleted copy constructor.
     */
    KineticMcChainAbstract(const KineticMcChainAbstract &) = delete;

    /**
     * @brief Deleted assignment operator.
     */
    void operator=(const mc::KineticMcChainAbstract &) = delete;

  protected:
    /**
     * @brief Performs one step of the simulation.
     */
    void OneStepSimulation() override;

    // helpful properties

    /**
     * @brief Lattice ID of the previous jump site.
     *
     */
    size_t previous_j_lattice_id_;

    /**
     * @brief Total rate for jump to nearest neighbours.
     *
     */
    double total_rate_i_{0.0};

    /**
     * @brief Neighbours of choosen lattice Id for the jump.
     *        j -> k -> i ->l
     *        k : current position
     *
     */
    vector<size_t> l_lattice_id_list_{};

    MPI_Op mpi_op_{};
    MPI_Datatype mpi_datatype_{};
  };

  /**
   * @struct MpiData
   * @brief Structure for storing data used in MPI custom reduction.
   */
  struct MpiData
  {
    double beta_bar_k{0.0};
    double beta_k{0.0};
    double gamma_bar_k_j{0.0};
    double gamma_k_j{0.0};
    double beta_k_j{0.0};
    double alpha_k_j{0.0};
    double ts_numerator{0.0};
    double ts_j_numerator{0.0};
  };

  /**
   * @brief Performs a custom reduction operation for MpiData.
   *
   * @param input_buffer Input buffer.
   * @param output_buffer Output buffer.
   * @param len Number of elements.
   * @param datatype MPI datatype.
   */
  inline void DataSum(void *input_buffer,
                      void *output_buffer,
                      int *len,
                      [[maybe_unused]] MPI_Datatype *datatype)
  {
    auto *input = static_cast<MpiData *>(input_buffer);
    auto *output = static_cast<MpiData *>(output_buffer);
    for (int i = 0; i < *len; ++i)
    {
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

  inline void DefineStruct(MPI_Datatype *datatype)
  {
    const int count = 8;
    int block_lens[count];
    MPI_Datatype types[count];
    MPI_Aint displacements[count];

    for (int i = 0; i < count; i++)
    {
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

#endif // LMC_LMC_MC_INCLUDE_KINETICMCABSTRACT_H_
