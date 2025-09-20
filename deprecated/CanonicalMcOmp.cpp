/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! \file CanonicalMcOmp.cpp
 *  @brief File for CanonicalMcOmp class implementation.
 */

#include "CanonicalMcOmp.h"
namespace mc
{
  CanonicalMcOmp::CanonicalMcOmp(Config config,
                                 unsigned long long int logDumpSteps,
                                 unsigned long long int configDumpStep,
                                 unsigned long long int maximumSteps,
                                 unsigned long long int thermodynamicAveragingSteps,
                                 unsigned long long int restartSteps,
                                 double restartEnergy,
                                 double temperature,
                                 EnergyPredictor &energyChangePredictor)
      : CanonicalMcAbstract(move(config),
                            logDumpSteps,
                            configDumpStep,
                            maximumSteps,
                            thermodynamicAveragingSteps,
                            restartSteps,
                            restartEnergy,
                            temperature,
                            energyChangePredictor)
  {
    if (world_size_ != 1)
    {
      cout << "Must use 1 process. Terminating...\n"
           << endl;
      MPI_Finalize();
      exit(0);
    }
#pragma omp parallel default(none) shared(cout)
    {
#pragma omp master
      {
        cout << "Using " << omp_get_num_threads() << " threads." << endl;
        num_threads_ = static_cast<size_t>(omp_get_num_threads());
      }
    }
    event_vector_.reserve(static_cast<size_t>(omp_get_num_threads()));
  }
  void CanonicalMcOmp::BuildEventVector()
  {
    event_vector_.clear();
    unavailable_position_.clear();

    pair<size_t, size_t> lattice_id_jump_pair;
    for (size_t i = 0; i < num_threads_; ++i)
    {
      size_t ct = 0;
      do
      {
        lattice_id_jump_pair = GenerateLatticeSiteIdJumpPair();
        ct++;
        if (ct == 50)
        {
          break;
        }
      } while (unavailable_position_.find(lattice_id_jump_pair.first) != unavailable_position_.end() || unavailable_position_.find(lattice_id_jump_pair.second) != unavailable_position_.end());
      if (ct == 50)
      {
        break;
      }
      for (auto selected_lattice_index : {lattice_id_jump_pair.first, lattice_id_jump_pair.second})
      {
        unavailable_position_.emplace(selected_lattice_index);

        auto firstNN = config_.GetNeighborLatticeIdVectorOfLattice(selected_lattice_index, 1);
        copy(firstNN.begin(), firstNN.end(),
             inserter(unavailable_position_,
                      unavailable_position_.end()));

        auto secondNN = config_.GetNeighborLatticeIdVectorOfLattice(selected_lattice_index, 2);
        copy(secondNN.begin(), secondNN.end(),
             inserter(unavailable_position_,
                      unavailable_position_.end()));

        auto thirdNN = config_.GetNeighborLatticeIdVectorOfLattice(selected_lattice_index, 3);
        copy(thirdNN.begin(), thirdNN.end(),
             inserter(unavailable_position_,
                      unavailable_position_.end()));
      }
      event_vector_.emplace_back(lattice_id_jump_pair, 0);
    }
#pragma omp parallel for default(none)
    for (auto &event : event_vector_)
    {
      event.second = energyChangePredictor_.GetEnergyChange(config_, event.first);
    }
  }
  void CanonicalMcOmp::Simulate()
  {

    while (steps_ <= maximumSteps_)
    {
      BuildEventVector();
      for (auto [lattice_id_jump_pair, dE] : event_vector_)
      {
        thermodynamicAveraging_.AddEnergy(energy_);
        Dump();
        SelectEvent(lattice_id_jump_pair, dE);
        ++steps_;
      }
    }
  }
} // mc