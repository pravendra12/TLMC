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
  CanonicalMcOmp::CanonicalMcOmp(TiledSupercell tiledSupercell,
                                 unsigned long long int logDumpSteps,
                                 unsigned long long int configDumpStep,
                                 unsigned long long int maximumSteps,
                                 unsigned long long int thermodynamicAveragingSteps,
                                 unsigned long long int restartSteps,
                                 double restartEnergy,
                                 double temperature,
                                 EnergyPredictorTLMC &energyChangePredictor)
      : CanonicalMcAbstract(move(tiledSupercell),
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

    pair<LatticeSiteMapping, LatticeSiteMapping> latticeIdJumpPair;
    size_t cubeIdx;
    for (size_t i = 0; i < num_threads_; ++i)
    {
      size_t ct = 0;
      do
      {
        // lattice_id_jump_pair = GenerateLatticeSiteIdJumpPair();
        latticeIdJumpPair = GenerateJumpPairInCube();
        cubeIdx = latticeIdJumpPair.first.smallConfigId;

        ct++;
        if (ct == 50)
        {
          break;
        }
      } while (unavailable_position_.find(cubeIdx) != unavailable_position_.end());
      if (ct == 50)
      {
        break;
      }

      unavailable_position_.emplace(cubeIdx);

      // const std::vector<size_t> &Cube::GetNeighbors(size_t siteIndex)
      auto neighboursOfCube = tiledSupercell_.GetCube().GetNeighbors(cubeIdx);

      copy(neighboursOfCube.begin(), neighboursOfCube.end(),
           inserter(unavailable_position_,
                    unavailable_position_.end()));

      event_vector_.emplace_back(latticeIdJumpPair, 0);
    }
#pragma omp parallel for default(none)
    for (auto &event : event_vector_)
    {
      event.second = energyChangePredictor_.GetEnergyChange(tiledSupercell_, event.first);
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