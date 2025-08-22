#include "SimulatedAnnealing.h"
namespace mc
{
  SimulatedAnnealing::SimulatedAnnealing(Config config,
                                         Config supercellConfig,
                                         const unsigned long long int logDumpSteps,
                                         const unsigned long long int configDumpSteps,
                                         const unsigned long long int maximumSteps,
                                         const unsigned long long int restartSteps,
                                         const double restartEnergy,
                                         const double initialTemperature,
                                         const set<Element> &elementSet,
                                         const string &predictorFilename)
      : McAbstract(move(config),
                   supercellConfig,
                   logDumpSteps,
                   configDumpSteps,
                   maximumSteps,
                   0, // thermodynamic_averaging_steps
                   restartSteps,
                   restartEnergy,
                   0,                  // restart time
                   initialTemperature, // temperature
                   elementSet,
                   predictorFilename,
                   "sa_log.txt"),
        energy_change_predictor_(predictorFilename,
                                 config,
                                 supercellConfig,
                                 elementSet),
        atom_index_selector_(0, config_.GetNumAtoms() - 1),
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
      // config_.WriteConfig(to_string(steps_) + ".cfg", config_);
      config_.WriteConfig(to_string(steps_) + ".cfg.gz", config_);
    }
    if (steps_ == maximumSteps_)
    {
      config_.WriteConfig("end.cfg", config_);
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
           << energy_ / double(config_.GetNumAtoms()) << endl;
    }
  }

  void SimulatedAnnealing::UpdateTemperature()
  {
    static const double alpha = 1. - 3. / static_cast<double>(maximumSteps_);
    temperature_ = initialTemperature_ * pow(alpha, steps_);
    beta_ = 1.0 / constants::kBoltzmann / temperature_;
  }

  pair<size_t, size_t> SimulatedAnnealing::GenerateLatticeIdJumpPair()
  {
    size_t lattice_id1, lattice_id2;
    do
    {
      lattice_id1 = atom_index_selector_(generator_);
      lattice_id2 = atom_index_selector_(generator_);
    } while (config_.GetElementOfLattice(lattice_id1) == config_.GetElementOfLattice(lattice_id2));
    return {lattice_id1, lattice_id2};
  }

  void SimulatedAnnealing::SelectEvent(const pair<size_t, size_t> &lattice_id_jump_pair,
                                       const double dE)
  {
    if (dE < 0)
    {
      config_.LatticeJump(lattice_id_jump_pair);
      energy_ += dE;
      absolute_energy_ += dE;
    }
    else
    {
      double possibility = exp(-dE * beta_);
      double random_number = unitDistribution_(generator_);
      if (random_number < possibility)
      {
        config_.LatticeJump(lattice_id_jump_pair);
        energy_ += dE;
        absolute_energy_ += dE;
      }
    }
  }

  void SimulatedAnnealing::Simulate()
  {
    while (steps_ <= maximumSteps_)
    {
      auto lattice_id_jump_pair = GenerateLatticeIdJumpPair();
      auto dE = energy_change_predictor_.GetDeSwap(config_, lattice_id_jump_pair);
      Dump();
      UpdateTemperature();
      SelectEvent(lattice_id_jump_pair, dE);
      ++steps_;
    }
  }

}