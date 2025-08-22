/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/**
 * @file KineticMcAbstract.h
 * @brief File contains implementation of KineticMcAbstract Class.
 */

#include "KineticMcAbstract.h"

namespace mc
{

  KineticMcFirstAbstract::KineticMcFirstAbstract(Config config,
                                                 Config supercellConfig,
                                                 const unsigned long long int logDumpSteps,
                                                 const unsigned long long int configDumpSteps,
                                                 const unsigned long long int maximumSteps,
                                                 const unsigned long long int thermodynamicAveragingSteps,
                                                 const unsigned long long int restartSteps,
                                                 const double restartEnergy,
                                                 const double restartTime,
                                                 const double temperature,
                                                 const set<Element> &elementSet,
                                                 const string &predictorFilename,
                                                 const string &timeTemperatureFilename,
                                                 const bool isRateCorrector,
                                                 const Eigen::RowVector3d &vacancyTrajectory)
      : McAbstract(move(config),
                   supercellConfig,
                   logDumpSteps,
                   configDumpSteps,
                   maximumSteps,
                   thermodynamicAveragingSteps,
                   restartSteps,
                   restartEnergy,
                   restartTime,
                   temperature,
                   elementSet,
                   predictorFilename,
                   "kmc_log.txt"),
        kEventListSize_(config.GetNeighborLatticeIdVectorOfLattice(0, 1).size()),
        ceParams_(predictorFilename),
        vacancyMigrationPredictor_(ceParams_,
                                   config),
        energyChangePredictor_(predictorFilename,
                               config,
                               supercellConfig,
                               elementSet),
        timeTemperatureInterpolator_(timeTemperatureFilename),
        isTimeTemperatureInterpolator_(!timeTemperatureFilename.empty()),
        isRateCorrector_(isRateCorrector),
        vacancyLatticeId_(config_.GetVacancyLatticeId()),
        vacancyTrajectory_(vacancyTrajectory),
        event_k_i_list_(kEventListSize_)
  {
  }

  KineticMcFirstAbstract::~KineticMcFirstAbstract() = default;

  void KineticMcFirstAbstract::UpdateTemperature()
  {
    if (isTimeTemperatureInterpolator_)
    {
      temperature_ = timeTemperatureInterpolator_.GetTemperature(time_);
      beta_ = 1.0 / constants::kBoltzmann / temperature_;
    }
  }

  //
  // double KineticMcFirstAbstract::GetTimeCorrectionFactor() {
  //   if (isRateCorrector_) {
  //     return rate_corrector_.GetTimeCorrectionFactor(temperature_);
  //   }
  //   return 1.0;
  // }

  void KineticMcFirstAbstract::Dump() const
  {
    if (is_restarted_)
    {
      is_restarted_ = false;
      return;
    }
    if (world_rank_ != 0)
    {
      return;
    }
    if (steps_ == 0)
    {
      ofs_ << "steps\ttime\taverage_time\ttemperature\tenergy\taverage_energy\tEa\tdE\tEa_Backward\tselected\tvacancy_trajectory";
      ofs_ << endl;
    }
    if (steps_ % configDumpSteps_ == 0)
    {
      config_.WriteConfig(to_string(steps_) + ".cfg.gz", config_);
    }
    if (steps_ == maximumSteps_)
    {
      config_.WriteConfig("end.cfg", config_);
    }

    unsigned long long int logDumpSteps;
    if (steps_ > 10 * logDumpSteps_)
    {
      logDumpSteps = logDumpSteps_;
    }
    else
    {
      logDumpSteps = static_cast<unsigned long long int>(
          pow(10, static_cast<unsigned long long int>(log10(steps_ + 1) - 1)));
      logDumpSteps = max(logDumpSteps, static_cast<unsigned long long int>(1));
      logDumpSteps = min(logDumpSteps, logDumpSteps_);
    }
    if (steps_ % logDumpSteps == 0)
    {
      ofs_ << steps_ << '\t'
           << time_ << '\t'
           << 1 / (total_rate_k_ * constants::kPrefactor) << "\t"
           << temperature_ << '\t'
           << energy_ << '\t'
           << thermodynamicAveraging_.GetThermodynamicAverage(beta_) << "\t"
           << event_k_i_.GetForwardBarrier() << '\t'
           << event_k_i_.GetEnergyChange() << '\t'
           << event_k_i_.GetBackwardBarrier() << "\t"
           << event_k_i_.GetIdJumpPair().second << '\t'
           << vacancyTrajectory_ << endl;
    }
  }

  size_t KineticMcFirstAbstract::SelectEvent() const
  {
    const double random_number = unitDistribution_(generator_);
    auto it = lower_bound(
        event_k_i_list_.begin(), event_k_i_list_.end(), random_number, [](const auto &lhs, double value)
        { return lhs.GetCumulativeProbability() < value; });

    // If not find (maybe generated 1), which rarely happens, returns the last event
    if (it == event_k_i_list_.cend())
    {
      it--;
    }
    if (world_size_ > 1)
    {
      int event_id = static_cast<int>(it - event_k_i_list_.cbegin());
      MPI_Bcast(&event_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
      return static_cast<size_t>(event_id);
    }
    else
    {
      return static_cast<size_t>(it - event_k_i_list_.cbegin());
    }
  }

  void KineticMcFirstAbstract::Debug(double one_step_time) const
  {
    if (world_rank_ == 0)
    {
      if (isnan(one_step_time) or isinf(one_step_time) or one_step_time < 0.0)
      {

        config_.WriteConfig("debug" + to_string(steps_) + ".cfg.gz", config_);
        cerr << "Invalid time step: " << one_step_time << endl;
        cerr << "For each event: Energy Barrier, Energy Change, Probability, " << endl;
        for (auto &event : event_k_i_list_)
        {
          cerr << event.GetForwardBarrier() << '\t' << event.GetEnergyChange() << '\t'
               << event.GetCumulativeProbability() << endl;
        }
        throw runtime_error("Invalid time step");
      }
    }
  }

  void KineticMcFirstAbstract::OneStepSimulation()
  {
    UpdateTemperature();

    thermodynamicAveraging_.AddEnergy(energy_);

    BuildEventList();

    // double one_step_time = CalculateTime() * GetTimeCorrectionFactor();
    double one_step_time = CalculateTime();
    Debug(one_step_time);

    event_k_i_ = event_k_i_list_[SelectEvent()];

    Dump();

    // modify
    time_ += one_step_time;
    energy_ += event_k_i_.GetEnergyChange();
    absolute_energy_ += event_k_i_.GetEnergyChange();

    Eigen::RowVector3d vacancyTrajectory_step =
        (config_.GetRelativeDistanceVectorLattice(vacancyLatticeId_,
                                                  event_k_i_.GetIdJumpPair().second))
            .transpose() *
        config_.GetBasis();

    vacancyTrajectory_ += vacancyTrajectory_step;

    config_.LatticeJump(event_k_i_.GetIdJumpPair());
    ++steps_;
    vacancyLatticeId_ = event_k_i_.GetIdJumpPair().second;
  }

  void KineticMcFirstAbstract::Simulate()
  {
    while (steps_ <= maximumSteps_)
    {
      OneStepSimulation();
    }
  }

  KineticMcChainAbstract::KineticMcChainAbstract(Config config,
                                                 Config supercellConfig,
                                                 const unsigned long long int logDumpSteps,
                                                 const unsigned long long int configDumpSteps,
                                                 const unsigned long long int maximumSteps,
                                                 const unsigned long long int thermodynamicAveragingSteps,
                                                 const unsigned long long int restartSteps,
                                                 const double restartEnergy,
                                                 const double restartTime,
                                                 const double temperature,
                                                 const set<Element> &elementSet,
                                                 const string &predictorFilename,
                                                 const string &timeTemperatureFilename,
                                                 const bool isRateCorrector,
                                                 const Eigen::RowVector3d &vacancyTrajectory)
      : KineticMcFirstAbstract(move(config),
                               supercellConfig,
                               logDumpSteps,
                               configDumpSteps,
                               maximumSteps,
                               thermodynamicAveragingSteps,
                               restartSteps,
                               restartEnergy,
                               restartTime,
                               temperature,
                               elementSet,
                               predictorFilename,
                               timeTemperatureFilename,
                               isRateCorrector,
                               vacancyTrajectory),
        previous_j_lattice_id_(config_.GetNeighborLatticeIdVectorOfLattice(vacancyLatticeId_, 1)[0]),
        l_lattice_id_list_(kEventListSize_)
  {
    MPI_Op_create(DataSum, 1, &mpi_op_);
    DefineStruct(&mpi_datatype_);
  }

  void KineticMcChainAbstract::OneStepSimulation()
  {
    KineticMcFirstAbstract::OneStepSimulation();
    previous_j_lattice_id_ = event_k_i_.GetIdJumpPair().first;
  }

  KineticMcChainAbstract::~KineticMcChainAbstract()
  {
    MPI_Op_free(&mpi_op_);
    MPI_Type_free(&mpi_datatype_);
  }

} // namespace mc
