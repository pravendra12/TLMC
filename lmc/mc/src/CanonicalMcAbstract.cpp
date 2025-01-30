/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 9/26/23 8:10 PM                                                           *
 **************************************************************************************************/

#include "CanonicalMcAbstract.h"
#include <utility>
#include <chrono>
#include <omp.h>
namespace mc {
// static std::set<Element> GetElementSetFromSolventAndSolute(
//     Element solvent_element, const std::set<Element> &solute_element_set) {
//
//   std::set<Element> element_set;
//   element_set.insert(solvent_element);
//   for (const auto &solute_element: solute_element_set) {
//     element_set.insert(solute_element);
//   }
//   return element_set;
// }
CanonicalMcAbstract::CanonicalMcAbstract(Config config,
                                         Config supercell_config,
                                         const unsigned long long int log_dump_steps,
                                         const unsigned long long int config_dump_steps,
                                         const unsigned long long int maximum_steps,
                                         const unsigned long long int thermodynamic_averaging_steps,
                                         const unsigned long long int restart_steps,
                                         const double restart_energy,
                                         const double temperature,
    // const double initial_temperature,
    // const double decrement_temperature,
                                         const std::set<Element> &element_set,
                                         const size_t max_cluster_size,
                                         const size_t max_bond_order,
                                         const std::string &json_coefficients_filename)
    : McAbstract(std::move(config),
                 supercell_config,
                 log_dump_steps,
                 config_dump_steps,
                 maximum_steps,
                 thermodynamic_averaging_steps,
                 restart_steps,
                 restart_energy,
                 0,
                 temperature,
                 element_set,
                 max_cluster_size,
                 max_bond_order,
                 json_coefficients_filename,
                 "cmc_log.txt"),
    // initial_temperature_(initial_temperature),
    // decrement_temperature_(decrement_temperature),
      energy_change_predictor_(json_coefficients_filename,
                               config,
                               supercell_config,
                               element_set,
                               max_cluster_size,
                               max_bond_order),
      vacancy_migration_predictor_(json_coefficients_filename),
      atom_index_selector_(0, config_.GetNumAtoms() - 1){
}

std::pair<size_t, size_t> CanonicalMcAbstract::GenerateLatticeIdJumpPair() {
  size_t lattice_id1, lattice_id2;
  do {
    lattice_id1 = atom_index_selector_(generator_);
    lattice_id2 = atom_index_selector_(generator_);
  } while (config_.GetElementOfLattice(lattice_id1)
      == config_.GetElementOfLattice(lattice_id2));
  return {lattice_id1, lattice_id2};
}

std::pair<size_t, size_t> CanonicalMcAbstract::GenerateVacancyAtomJumpPair() {
  size_t vacancy_id, lattice_id;
  do {
    vacancy_id = config_.GetVacancyLatticeId();
    lattice_id = atom_index_selector_(generator_);
  } while (config_.GetElementOfLattice(vacancy_id)
      == config_.GetElementOfLattice(lattice_id));
  return {vacancy_id, lattice_id};
}



// void CanonicalMcAbstract::UpdateTemperature() {
//   if (steps_ % maximum_steps_ == 0 && steps_ != 0) {
//     config_.WriteConfig("end_" + std::to_string(static_cast<int>(temperature_)) + "K.cfg",
//                         false);
//     temperature_ -= decrement_temperature_;
//     beta_ = 1.0 / constants::kBoltzmann / temperature_;
//   }
// }

void CanonicalMcAbstract::Dump() const {
  if (is_restarted_) {
    is_restarted_ = false;
    return;
  }
  if (steps_ == 0) {
    // config_.WriteLattice("lattice.txt");
    // config_.WriteElement("element.txt");
    // config_.WriteConfig("start.cfg",config_);
    // ofs_ << "steps\ttemperature\tenergy\taverage_energy\tabsolute_energy" << std::endl;
    ofs_ << "steps\ttemperature\tenergy\tenergy_per_atom" << std::endl;

  }
  if (steps_ % config_dump_steps_ == 0) {
    //config_.WriteMap("map" + std::to_string(steps_) + ".txt");
    config_.WriteConfig(std::to_string(steps_) + ".cfg", config_);
  }
  if (steps_ == maximum_steps_) {
    config_.WriteConfig("end.cfg", config_);
  }
  unsigned long long int log_dump_steps;
  if (steps_ > 10 * log_dump_steps_) {
    log_dump_steps = log_dump_steps_;
  } else {
    log_dump_steps = static_cast<unsigned long long int>(
        std::pow(10, static_cast<unsigned long long int>(std::log10(
            steps_ + 1) - 1)));
    log_dump_steps = std::max(log_dump_steps, static_cast<unsigned long long int>(1));
    log_dump_steps = std::min(log_dump_steps, log_dump_steps_);
  }
  if (steps_ % log_dump_steps == 0) {
    ofs_ << steps_ << '\t' << temperature_ << '\t' << energy_ << '\t'
         << energy_/config_.GetNumAtoms() << std::endl;
  }
}
void CanonicalMcAbstract::SelectEvent(const std::pair<size_t, size_t> &lattice_id_jump_pair,
                                      const double dE) {
  if (dE < 0) {

    // std::cout << "dE is less than 0" << std::endl;

    config_.LatticeJump(lattice_id_jump_pair);
    energy_ += dE;
    absolute_energy_ += dE;
    
  } else {
    double possibility = std::exp(-dE * beta_);
    double random_number = unit_distribution_(generator_);
    if (random_number < possibility) {
      config_.LatticeJump(lattice_id_jump_pair);
      energy_ += dE;
      absolute_energy_ += dE;
    }
  }

  // std::cout << "Absolute Energy: " << absolute_energy_ << std::endl;

}
} // mc