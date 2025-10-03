/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! \file CanonicalMcAbstract.cpp
 *  @brief File for CanonicalMcAbstract class implementation.
 */

#include "CanonicalMcAbstract.h"

namespace mc
{
  CanonicalMcAbstract::CanonicalMcAbstract(TiledSupercell tiledSupercell,
                                           const unsigned long long int logDumpSteps,
                                           const unsigned long long int configDumpSteps,
                                           const unsigned long long int maximumSteps,
                                           const unsigned long long int thermodynamicAveragingSteps,
                                           const unsigned long long int restartSteps,
                                           const double restartEnergy,
                                           const double temperature,
                                           EnergyPredictorTLMC &energyChangePredictor)
      : McAbstract(move(tiledSupercell),
                   logDumpSteps,
                   configDumpSteps,
                   maximumSteps,
                   thermodynamicAveragingSteps,
                   restartSteps,
                   restartEnergy,
                   0,
                   temperature,
                   "cmc_log.txt"),
        energyChangePredictor_(energyChangePredictor),
        atomIndexSelector_(0, tiledSupercell_.GetTotalNumOfSites() - 1),
        cubeIndexSelector_(0, tiledSupercell_.GetNumOfSmallConfig() - 1),
        smallConfigAtomIndexSelector_(0, tiledSupercell.GetNumOfSitesPerSmallConfig()-1)
  {
  }

  pair<LatticeSiteMapping, LatticeSiteMapping> CanonicalMcAbstract::GenerateLatticeSiteIdJumpPair()
  {
    LatticeSiteMapping latticeSiteId1;
    LatticeSiteMapping latticeSiteId2;
    do
    {
      auto atomId1 = atomIndexSelector_(generator_);
      latticeSiteId1 = tiledSupercell_.GetLatticeSiteMappingFromAtomId(atomId1);

      auto atomId2 = atomIndexSelector_(generator_);
      latticeSiteId2 = tiledSupercell_.GetLatticeSiteMappingFromAtomId(atomId2);

    } while (tiledSupercell_.GetElementAtSite(latticeSiteId1) == tiledSupercell_.GetElementAtSite(latticeSiteId2));

    return {latticeSiteId1, latticeSiteId2};
  }

  pair<LatticeSiteMapping, LatticeSiteMapping> CanonicalMcAbstract::GenerateJumpPairInCube()
  {
    auto randomCubeIdx = cubeIndexSelector_(generator_);

    LatticeSiteMapping latticeSiteId1;
    LatticeSiteMapping latticeSiteId2;
    do
    {
      auto atomId1 = smallConfigAtomIndexSelector_(generator_);
      latticeSiteId1 = LatticeSiteMapping(atomId1, randomCubeIdx);

      auto atomId2 = smallConfigAtomIndexSelector_(generator_);
      latticeSiteId2 = LatticeSiteMapping(atomId2, randomCubeIdx);

    } while (tiledSupercell_.GetElementAtSite(latticeSiteId1) == tiledSupercell_.GetElementAtSite(latticeSiteId2));

    return {latticeSiteId1, latticeSiteId2};
  }

  pair<LatticeSiteMapping, LatticeSiteMapping> CanonicalMcAbstract::GenerateVacancyLatticeSiteIdJumpPair()
  {
    auto latticeSiteId1 = tiledSupercell_.GetVacancySiteId();

    LatticeSiteMapping latticeSiteId2;
    do
    {
      auto atomId2 = atomIndexSelector_(generator_);
      latticeSiteId2 = tiledSupercell_.GetLatticeSiteMappingFromAtomId(atomId2);
    } while (tiledSupercell_.GetElementAtSite(latticeSiteId1) == tiledSupercell_.GetElementAtSite(latticeSiteId2));

    return {latticeSiteId1, latticeSiteId2};
  }

  void CanonicalMcAbstract::Dump() const
  {
    if (is_restarted_)
    {
      is_restarted_ = false;
      return;
    }
    if (steps_ == 0)
    {
      ofs_ << "steps\ttemperature\tenergy\tenergy_per_atom\taverage_energy" << endl;
    }
    if (steps_ % configDumpSteps_ == 0)
    {
      // tiledSupercell_.WriteAtomVectorInfoToFile(to_string(steps_) + ".txt");
      tiledSupercell_.WriteAtomicIndicesToFile(to_string(steps_) + ".bin.gz");
    }
    if (steps_ == maximumSteps_)
    {
      // tiledSupercell_.WriteAtomVectorInfoToFile(to_string(steps_) + ".txt");
      tiledSupercell_.WriteAtomicIndicesToFile(to_string(steps_) + ".bin.gz");
    }
    unsigned long long int logDumpSteps;
    if (steps_ > 10 * logDumpSteps_)
    {
      logDumpSteps = logDumpSteps_;
    }
    else
    {
      logDumpSteps = static_cast<unsigned long long int>(
          pow(10, static_cast<unsigned long long int>(log10(
                                                          steps_ + 1) -
                                                      1)));
      logDumpSteps = max(logDumpSteps, static_cast<unsigned long long int>(1));
      logDumpSteps = min(logDumpSteps, logDumpSteps_);
    }
    if (steps_ % logDumpSteps == 0)
    {
      ofs_ << steps_ << '\t'
           << temperature_ << '\t'
           << energy_ << '\t'
           << energy_ / static_cast<double>(tiledSupercell_.GetTotalNumOfSites()) << "\t"
           << thermodynamicAveraging_.GetThermodynamicAverage(beta_)
           << endl;
    }
  }
  void CanonicalMcAbstract::SelectEvent(const pair<LatticeSiteMapping, LatticeSiteMapping> &lattice_id_jump_pair,
                                        const double dE)
  {
    if (dE < 0)
    {
      tiledSupercell_.LatticeJump(lattice_id_jump_pair);
      energy_ += dE;
      absolute_energy_ += dE;

      thermodynamicAveraging_.AddEnergy(dE);
    }
    else
    {
      double possibility = exp(-dE * beta_);
      double random_number = unitDistribution_(generator_);
      if (random_number < possibility)
      {
        tiledSupercell_.LatticeJump(lattice_id_jump_pair);
        energy_ += dE;
        absolute_energy_ += dE;

        thermodynamicAveraging_.AddEnergy(dE);
      }
    }
  }
} // mc