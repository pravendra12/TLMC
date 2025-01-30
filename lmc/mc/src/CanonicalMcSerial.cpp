/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 10/30/23 3:09 PM                                                          *
 **************************************************************************************************/

#include "CanonicalMcSerial.h"
#include <utility>
#include <chrono>
#include <omp.h>
namespace mc {
CanonicalMcSerial::CanonicalMcSerial(Config config,
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
    : CanonicalMcAbstract(std::move(config),
                          supercell_config,
                          log_dump_steps,
                          config_dump_steps,
                          maximum_steps,
                          thermodynamic_averaging_steps,
                          restart_steps,
                          restart_energy,
                          temperature,
    // initial_temperature,
    // decrement_temperature,
                          element_set,
                          max_cluster_size,
                          max_bond_order,
                          json_coefficients_filename) {
  if (world_size_ != 1) {
    std::cout << "Must use 1 precesses. Terminating...\n" << std::endl;
    MPI_Finalize();
    exit(0);
  }
#pragma omp parallel  default(none) shared(std::cout)
  {
#pragma omp master
    {
      std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
    }
  }
}

void CanonicalMcSerial::Simulate() {
  // std::cout << "Absolute Energy: " << absolute_energy_ << std::endl;
  // while (steps_ <= maximum_steps_ * static_cast<unsigned long long int>(initial_temperature_ / decrement_temperature_ + 1)) {
  while (steps_ <= maximum_steps_) {

    // Generate the lattice id jump pair
    // auto lattice_id_jump_pair = GenerateLatticeIdJumpPair();

    // Generate the vacancy and lattice id jump pair
    // Used for xSv determinant
    auto vacancy_lattice_id_jump_pair = GenerateVacancyAtomJumpPair();
    
    // This need to be derived using the  cluster expansion hamiltonian
    // As xSv method specifically predicts the dE lattice id jump pair which 
    // involves a vacancy, may be for binary alloy its not reliable as it only     
    // consider the first NN
    
    // For CE method
    // auto dE = energy_change_predictor_.GetDe(config_, lattice_id_jump_pair);

    // Using xSv Determinant
    auto dE = vacancy_migration_predictor_.GetDiff(config_, 
                                                   vacancy_lattice_id_jump_pair);

    // std::cout << "{" << vacancy_lattice_id_jump_pair.first << "\t" 
    //                  << vacancy_lattice_id_jump_pair.second << "}\t" 
    //           << dE << std::endl;

    Dump();

    // SelectEvent(lattice_id_jump_pair, dE);
    SelectEvent(vacancy_lattice_id_jump_pair, dE);

    ++steps_;
  }
}
} // mc
