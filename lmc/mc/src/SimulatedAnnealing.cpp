#include "SimulatedAnnealing.h"

SimulatedAnnealing::SimulatedAnnealing(Config config,
                                       Config supercell_config,
                                       const unsigned long long int log_dump_steps,
                                       const unsigned long long int config_dump_steps,
                                       const unsigned long long int maximum_steps,
                                       const unsigned long long int restart_steps,
                                       const double restart_energy,
                                       const double initial_temperature,
                                       const std::set<Element> &element_set,                                   
                                       const size_t max_cluster_size,
                                       const size_t max_bond_order,
                                       const std::string &json_coefficients_filename)
    :McAbstract(std::move(config),
                 supercell_config,
                 log_dump_steps,
                 config_dump_steps,
                 maximum_steps,
                 0, // thermodynamic_averaging_steps,
                 restart_steps,
                 restart_energy,
                 0,
                 initial_temperature, //temperature,
                 element_set,
                 max_cluster_size,
                 max_bond_order,
                 json_coefficients_filename,
                 "sa_log.txt"),
      energy_change_predictor_(json_coefficients_filename,
                               config,
                               supercell_config,
                               element_set,
                               max_cluster_size,
                               max_bond_order),
      atom_index_selector_(0, config_.GetNumAtoms() - 1),
      initial_temperature_(initial_temperature)

      
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
         << "energy_per_atom" << std::endl;
  }
  if (steps_ % config_dump_steps_ == 0) 
  {
    config_.WriteConfig(std::to_string(steps_) + ".cfg", config_);
  }
  if (steps_ == maximum_steps_) 
  {
    config_.WriteConfig("end.cfg", config_);
  }
  
  // log steps
  unsigned long long int log_dump_steps;
  if (steps_ > 10 * log_dump_steps_) 
  {
    log_dump_steps = log_dump_steps_;
  } 
  else 
  {
    log_dump_steps = static_cast<unsigned long long int>(
        std::pow(10, static_cast<unsigned long long int>(std::log10(
            steps_ + 1) - 1)));
    log_dump_steps = std::max(log_dump_steps, static_cast<unsigned long long int>(1));
    log_dump_steps = std::min(log_dump_steps, log_dump_steps_);
  }
  
  if (steps_ > 0 && steps_ % log_dump_steps == 0) 
  {
      ofs_ << steps_ << '\t' 
          << temperature_ << '\t' 
          << energy_ << '\t'
          << energy_ / double(config_.GetNumAtoms()) << std::endl;
  }
}

void SimulatedAnnealing::UpdateTemperature() 
{
  static const double alpha = 1. - 3. / static_cast<double>(maximum_steps_);
  temperature_ = initial_temperature_ * std::pow(alpha, steps_);
  beta_ = 1.0 / constants::kBoltzmann / temperature_;
}

std::pair<size_t, size_t> SimulatedAnnealing::GenerateLatticeIdJumpPair() 
{
  size_t lattice_id1, lattice_id2;
  do 
  {
    lattice_id1 = atom_index_selector_(generator_);
    lattice_id2 = atom_index_selector_(generator_);
  } 
  while (config_.GetElementOfLattice(lattice_id1)
      == config_.GetElementOfLattice(lattice_id2));
  return {lattice_id1, lattice_id2};
}

void SimulatedAnnealing::SelectEvent(const std::pair<size_t, size_t> &lattice_id_jump_pair,
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
    double possibility = std::exp(-dE * beta_);
    double random_number = unit_distribution_(generator_);
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
  while (steps_ <= maximum_steps_) {
    auto lattice_id_jump_pair = GenerateLatticeIdJumpPair();
    auto dE = energy_change_predictor_.GetDe(config_, lattice_id_jump_pair);
    Dump();
    UpdateTemperature();
    SelectEvent(lattice_id_jump_pair, dE);
    ++steps_;
  }
}


