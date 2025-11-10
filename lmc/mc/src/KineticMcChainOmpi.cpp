/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/**
 * @file KineticMcChainOmpi.h
 * @brief Implementation of the KineticMcChainOmpi class for OMPI-based kinetic Monte Carlo simulations.
 */

#include "KineticMcChainOmpi.h"

namespace mc
{

  KineticMcChainOmpi::KineticMcChainOmpi(TiledSupercell tiledSupercell,
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
                                         const Eigen::RowVector3d &vacancyTrajectory)
      : KineticMcChainAbstract(move(tiledSupercell),
                               logDumpSteps,
                               configDumpSteps,
                               maximumSteps,
                               thermodynamicAveragingSteps,
                               restartSteps,
                               restartEnergy,
                               restartTime,
                               temperature,
                               vacancyMigrationPredictor,
                               timeTemperatureFilename,
                               rateCorrector,
                               vacancyTrajectory)
  {
    if (world_size_ != kEventListSize_)
    {
      std::cout << "Must use " << kEventListSize_ << " processes. Terminating...\n"
                << std::endl;
      MPI_Finalize();
      exit(0);
    }
    if (world_rank_ == 0)
    {
      std::cout << "Using " << world_size_ << " processes." << std::endl;
#pragma omp parallel default(none) shared(std::cout)
      {
#pragma omp master
        {
          std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
        }
      }
    }
  }

  // Update  first_event_ki and l_index_list for each process
  void KineticMcChainOmpi::BuildEventList()
  {
    const auto k_lattice_id = vacancyLatticeSiteId_;

    // LATTICE_ID , ENCODED_CONFIG_IDX
    const auto encoded_i_lattice_site_id =
        tiledSupercell_.GetNeighborLatticeIdVectorOfLattice(k_lattice_id.latticeId, 1)[static_cast<size_t>(world_rank_)];

    const auto i_lattice_id = tiledSupercell_.GetLatticeSiteMappingFromEncoding(
        encoded_i_lattice_site_id,
        k_lattice_id);

    size_t it = 0;

    for (const auto &encoded_l_lattice_id : tiledSupercell_.GetNeighborLatticeIdVectorOfLattice(i_lattice_id.latticeId, 1))
    {
      const auto l_lattice_id = tiledSupercell_.GetLatticeSiteMappingFromEncoding(
          encoded_l_lattice_id,
          i_lattice_id);

      l_lattice_id_list_[it] = l_lattice_id;
      ++it;
    }

    total_rate_k_ = 0.0;
    total_rate_i_ = 0.0;

    tiledSupercell_.LatticeJump({k_lattice_id, i_lattice_id});

#pragma omp parallel default(none) shared(i_lattice_id, k_lattice_id) reduction(+ : total_rate_i_)
    {
#pragma omp for
      for (size_t ii = 0; ii < kEventListSize_; ++ii)
      {

        const auto l_lattice_id = l_lattice_id_list_[ii];

        JumpEvent event_i_l(
            // Jump Pair
            // {Vacancy, Migrating Atom}
            {i_lattice_id, l_lattice_id},

            vacancyMigrationPredictor_.GetBarrierAndDeltaE(tiledSupercell_,
                                                           {i_lattice_id,
                                                            l_lattice_id}),

            // Thermodynamics Beta
            beta_);

        if (l_lattice_id == k_lattice_id)
        {
          event_k_i_ = event_i_l.GetReverseJumpEvent();
        }
        auto r_i_l = event_i_l.GetForwardRate();
        total_rate_i_ += r_i_l;
      }
    }
    tiledSupercell_.LatticeJump({i_lattice_id, k_lattice_id});

    double rate_k_i = event_k_i_.GetForwardRate();
    MPI_Allreduce(&rate_k_i, &total_rate_k_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    event_k_i_.CalculateProbability(total_rate_k_);
  }

  double KineticMcChainOmpi::CalculateTime()
  {
    const auto probability_k_i = event_k_i_.GetProbability();
    const auto probability_i_k = event_k_i_.GetBackwardRate() / total_rate_i_;
    // Probability that flicker event happens, vacancy jump backwards
    const double beta_bar_k_i = probability_k_i * probability_i_k;
    // Probability that flicker event does not happen, vacancy jump forwards
    const double beta_k_i = probability_k_i * (1 - probability_i_k);

    bool is_previous_event = event_k_i_.GetIdJumpPair().second == previous_j_lattice_id_;

    // Time in first order Kmc, same for all
    const double t_1 = 1 / total_rate_k_ / constants::kPrefactor;
    // Time in neighbors, different for each process
    const double t_i = 1 / total_rate_i_ / constants::kPrefactor;
    MpiData mpi_data_helper{
        beta_bar_k_i,
        beta_k_i,
        is_previous_event ? 0.0 : beta_bar_k_i,
        is_previous_event ? 0.0 : beta_k_i,
        is_previous_event ? beta_k_i : 0.0,
        is_previous_event ? probability_k_i : 0.0,
        (t_1 + t_i) * beta_bar_k_i,
        is_previous_event ? 0.0 : (t_1 + t_i) * beta_bar_k_i};
    MpiData mpi_data{};
    MPI_Allreduce(&mpi_data_helper, &mpi_data, 1, mpi_datatype_, mpi_op_, MPI_COMM_WORLD);

    const auto beta_bar_k = mpi_data.beta_bar_k;
    const auto beta_k = mpi_data.beta_k;
    // sum of beta_bar_k_i such that i != j
    const auto gamma_bar_k_j = mpi_data.gamma_bar_k_j;
    // sum of beta_k_i such that i = j
    const auto gamma_k_j = mpi_data.gamma_k_j;
    const auto beta_k_j = mpi_data.beta_k_j;
    const auto alpha_k_j = mpi_data.alpha_k_j;
    const auto ts = mpi_data.ts_numerator / beta_bar_k;
    const auto ts_j = mpi_data.ts_j_numerator / gamma_bar_k_j;

    const double one_over_one_minus_a_j = 1 / (1 - alpha_k_j);
    const double t_2 = one_over_one_minus_a_j * (gamma_k_j * t_1 + gamma_bar_k_j * (ts_j + t_1 + beta_bar_k / beta_k * ts));

    double second_order_probability;
    if (is_previous_event)
    {
      second_order_probability = one_over_one_minus_a_j * (gamma_bar_k_j / beta_k) * beta_k_j;
    }
    else
    {
      second_order_probability = one_over_one_minus_a_j * (1 + gamma_bar_k_j / beta_k) * beta_k_i;
    }
    event_k_i_.SetProbability(second_order_probability);
    MPI_Allgather(&event_k_i_,
                  sizeof(JumpEvent),
                  MPI_BYTE,
                  event_k_i_list_.data(),
                  sizeof(JumpEvent),
                  MPI_BYTE,
                  MPI_COMM_WORLD);
    double cumulative_probability = 0.0;
    for (auto &event : event_k_i_list_)
    {
      cumulative_probability += event.GetProbability();
      event.SetCumulativeProbability(cumulative_probability);
    }
    return t_2;
  }
} // mc