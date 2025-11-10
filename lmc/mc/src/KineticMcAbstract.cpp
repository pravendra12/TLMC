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

  KineticMcFirstAbstract::KineticMcFirstAbstract(TiledSupercell tiledSupercell,
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
      : McAbstract(move(tiledSupercell),
                   logDumpSteps,
                   configDumpSteps,
                   maximumSteps,
                   thermodynamicAveragingSteps,
                   restartSteps,
                   restartEnergy,
                   restartTime,
                   temperature,
                   "kmc_log.txt"),
        kEventListSize_(
            tiledSupercell_.GetNeighborLatticeIdVectorOfLattice(0, 1).size()),
        vacancyMigrationPredictor_(
            vacancyMigrationPredictor),
        timeTemperatureInterpolator_(
            timeTemperatureFilename),
        isTimeTemperatureInterpolator_(!timeTemperatureFilename.empty()),
        rateCorrector_(move(rateCorrector)),
        vacancyLatticeSiteId_(tiledSupercell_.GetVacancySiteId()),
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

  double KineticMcFirstAbstract::GetTimeCorrectionFactor()
  {
    if (rateCorrector_)
    {
      return rateCorrector_->GetTimeCorrectionFactor(temperature_);
    }
    return 1.0;
  }

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
      // config_.WriteConfig(to_string(steps_) + ".cfg.gz", config_);
      // tiledSupercell_.WriteAtomVectorInfoToFile(to_string(steps_) + ".txt");
      tiledSupercell_.WriteAtomicIndicesToFile(to_string(steps_) + ".bin.gz");
    }
    if (steps_ == maximumSteps_)
    {
      tiledSupercell_.WriteAtomicIndicesToFile(to_string(steps_) + ".bin.gz");
      // tiledSupercell_.WriteAtomVectorInfoToFile(to_string(steps_) + ".txt");
      // config_.WriteConfig("end.cfg", config_);
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
           << "(" << event_k_i_.GetIdJumpPair().second.latticeId
           << ", " << event_k_i_.GetIdJumpPair().second.smallConfigId << ")\t"
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

        // config_.WriteConfig("debug" + to_string(steps_) + ".cfg.gz", config_);
        tiledSupercell_.WriteAtomVectorInfoToFile("debug" + to_string(steps_) + ".txt");
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

    double one_step_time = CalculateTime() * GetTimeCorrectionFactor();
    // double one_step_time = CalculateTime();
    Debug(one_step_time);

    event_k_i_ = event_k_i_list_[SelectEvent()];

    Dump();

    // modify
    time_ += one_step_time;
    energy_ += event_k_i_.GetEnergyChange();
    absolute_energy_ += event_k_i_.GetEnergyChange();

    Eigen::RowVector3d vacancyTrajectory_step =
        (tiledSupercell_.GetRelativeDistanceVectorLattice(
             vacancyLatticeSiteId_,
             event_k_i_.GetIdJumpPair().second))
            .transpose() *
        tiledSupercell_.GetSuperBasis();

    vacancyTrajectory_ += vacancyTrajectory_step;

    tiledSupercell_.LatticeJump(event_k_i_.GetIdJumpPair());
    ++steps_;
    vacancyLatticeSiteId_ = event_k_i_.GetIdJumpPair().second;
  }

  void KineticMcFirstAbstract::Simulate()
  {
    while (steps_ <= maximumSteps_)
    {
      OneStepSimulation();
    }
  }

  KineticMcChainAbstract::KineticMcChainAbstract(TiledSupercell tiledSupercell,
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
      : KineticMcFirstAbstract(move(tiledSupercell),
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
                               vacancyTrajectory),
        previous_j_lattice_id_(
            tiledSupercell_.GetLatticeSiteMappingFromEncoding(
                tiledSupercell_.GetNeighborLatticeIdVectorOfLattice(vacancyLatticeSiteId_.latticeId, 1)[0],
                vacancyLatticeSiteId_)),
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
