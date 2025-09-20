#include "SimulatedAnnealing.h"
namespace mc
{
  SimulatedAnnealing::SimulatedAnnealing(TiledSupercell tiledSupercell,
                                         const unsigned long long int logDumpSteps,
                                         const unsigned long long int configDumpSteps,
                                         const unsigned long long int maximumSteps,
                                         const unsigned long long int restartSteps,
                                         const double restartEnergy,
                                         const double initialTemperature,
                                         EnergyPredictorTLMC &energyChangePredictor)
      : McAbstract(move(tiledSupercell),
                   logDumpSteps,
                   configDumpSteps,
                   maximumSteps,
                   0, // thermodynamic_averaging_steps
                   restartSteps,
                   restartEnergy,
                   0,                  // restart time
                   initialTemperature, // temperature
                   "sa_log.txt"),
        energyChangePredictor_(energyChangePredictor),
        atomIndexSelector_(0, tiledSupercell.GetTotalNumOfSites() - 1),
        initialTemperature_(initialTemperature)

  {
    ofs_.precision(16);
  }

  void SimulatedAnnealing::Dump() const
  {
    if (is_restarted_)
    {
      is_restarted_ = false;
      return;
    }
    // dump configuration
    if (steps_ == 0)
    {
      ofs_ << "steps\t"
           << "temperature\t"
           << "energy\t"
           << "energy_per_atom" << endl;
    }
    if (steps_ % configDumpSteps_ == 0)
    {
      tiledSupercell_.WriteAtomVectorInfoToFile(to_string(steps_) + ".txt");
    }
    if (steps_ == maximumSteps_)
    {
      tiledSupercell_.WriteAtomVectorInfoToFile(to_string(steps_) + ".txt");
    }

    // log steps
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

    if (steps_ > 0 && steps_ % logDumpSteps == 0)
    {
      ofs_ << steps_ << '\t'
           << temperature_ << '\t'
           << energy_ << '\t'
           << energy_ / double(tiledSupercell_.GetTotalNumOfSites()) << endl;
    }
  }

  void SimulatedAnnealing::UpdateTemperature()
  {
    static const double alpha = 1. - 3. / static_cast<double>(maximumSteps_);
    temperature_ = initialTemperature_ * pow(alpha, steps_);
    beta_ = 1.0 / constants::kBoltzmann / temperature_;
  }

  pair<LatticeSiteMapping, LatticeSiteMapping> SimulatedAnnealing::GenerateLatticeSiteIdJumpPair()
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

  pair<LatticeSiteMapping, LatticeSiteMapping> SimulatedAnnealing::GenerateVacancyLatticeSiteIdJumpPair()
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

  void SimulatedAnnealing::SelectEvent(const pair<LatticeSiteMapping, LatticeSiteMapping> &lattice_id_jump_pair,
                                       const double dE)
  {
    if (dE < 0)
    {
      tiledSupercell_.LatticeJump(lattice_id_jump_pair);
      energy_ += dE;
      absolute_energy_ += dE;
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
      }
    }
  }

  void SimulatedAnnealing::Simulate()
  {
    while (steps_ <= maximumSteps_)
    {
      // condition can provided based on the whether to use vacancy-lattice site pair
      // or lattice-lattice site pair
      auto lattice_id_jump_pair = GenerateLatticeSiteIdJumpPair();
      auto dE = energyChangePredictor_.GetEnergyChange(
          tiledSupercell_,
          lattice_id_jump_pair);
      Dump();
      UpdateTemperature();
      SelectEvent(lattice_id_jump_pair, dE);
      ++steps_;
    }
  }

}