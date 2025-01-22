/**************************************************************************************************
 * Copyright (c) 2023-2024. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 3/21/22 3:17 PM                                                                         *
 * @Last Modified by: pravendra12                                                                    *
 * @Last Modified time: 11/28/24 8:30 PM                                                            *
 **************************************************************************************************/

/*! \file  Home.h
 *  \brief File for the Home class implementation.
 */

#include "Home.h"

namespace api {

void Print(const Parameter &parameter) {
  std::cout << "Parameters" << std::endl;
  std::cout << "simulation_method: " << parameter.method << std::endl;

  if (parameter.method == "KineticMcFirstMpi" || parameter.method == "KineticMcFirstOmp" ||
      parameter.method == "KineticMcChainOmp" || parameter.method == "KineticMcChainMpi" ||
      parameter.method == "KineticMcChainOmpi") {

    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
              << std::endl;
    std::cout << "config_filename: " << parameter.config_filename_ << std::endl;
    std::cout << "map_filename: " << parameter.map_filename_ << std::endl;
    std::cout << "supercell_size: " << parameter.supercell_size_ << std::endl;
    std::cout << "lattice_param: " << parameter.lattice_param_ << std::endl;
    std::cout << "structure_type: " << parameter.structure_type_ << std::endl;
    std::cout << "log_dump_steps: " << parameter.log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << parameter.config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << parameter.maximum_steps_ << std::endl;
    std::cout << "max_cluster_size: " << parameter.max_cluster_size_ << std::endl;
    std::cout << "max_bond_order: " << parameter.max_bond_order_ << std::endl;
    std::cout << "cutoffs: ";
    std::transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                   std::ostream_iterator<std::string>(std::cout, " "),
                   [](auto cutoff) { return std::to_string(cutoff); });
    std::cout << std::endl;
    std::cout << "thermodynamic_averaging_steps: " << parameter.thermodynamic_averaging_steps_
              << std::endl;
    std::cout << "temperature: " << parameter.temperature_ << std::endl;
    std::cout << "restart_steps: " << parameter.restart_steps_ << std::endl;
    std::cout << "restart_energy: " << parameter.restart_energy_ << std::endl;
    std::cout << "restart_time: " << parameter.restart_time_ << std::endl;
    std::cout << "time_temperature_filename: " << parameter.time_temperature_filename_
              << std::endl;
    std::cout << "rate_corrector: " << parameter.rate_corrector_ << std::endl;
    std::cout << "vacancy_trajectory: " << parameter.vacancy_trajectory_
              << std::endl;
  } 
  
  else if (parameter.method == "SimulatedAnnealing") {
    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
              << std::endl;
    std::cout << "factor: " << parameter.factor_ << std::endl;
    std::cout << "solvent_element: " << parameter.solvent_element_ << std::endl;
    std::cout << "solute_element_set: ";
    std::copy(parameter.solute_element_set_.begin(), parameter.solute_element_set_.end(),
              std::ostream_iterator<std::string>(std::cout, " "));
    std::cout << std::endl;
    std::cout << "solute_number_set: ";
    std::transform(parameter.solute_number_set_.begin(), parameter.solute_number_set_.end(),
                   std::ostream_iterator<std::string>(std::cout, " "),
                   [](auto number) { return std::to_string(number); });
    std::cout << std::endl;
    std::cout << "log_dump_steps: " << parameter.log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << parameter.config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << parameter.maximum_steps_ << std::endl;
    std::cout << "early_stop_steps: " << parameter.early_stop_steps_ << std::endl;
    std::cout << "initial_temperature: " << parameter.initial_temperature_ << std::endl;
  } 
  
  else if (parameter.method == "CanonicalMcSerial" || parameter.method == "CanonicalMcOmp") {
    std::cout << "json_coefficients_filename: " << parameter.json_coefficients_filename_
              << std::endl;
    std::cout << "config_filename: " << parameter.config_filename_ << std::endl;
    std::cout << "map_filename: " << parameter.map_filename_ << std::endl;
    std::cout << "supercell_size: " << parameter.supercell_size_ << std::endl;
    std::cout << "lattice_param: " << parameter.lattice_param_ << std::endl;
    std::cout << "structure_type: " << parameter.structure_type_ << std::endl;
    std::cout << "max_cluster_size: " << parameter.max_cluster_size_ << std::endl;
    std::cout << "max_bond_order: " << parameter.max_bond_order_ << std::endl;
    std::cout << "cutoffs: ";
    std::transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                   std::ostream_iterator<std::string>(std::cout, " "),
                   [](auto cutoff) { return std::to_string(cutoff); });
    std::cout << std::endl;
    std::cout << "log_dump_steps: " << parameter.log_dump_steps_ << std::endl;
    std::cout << "config_dump_steps: " << parameter.config_dump_steps_ << std::endl;
    std::cout << "maximum_steps: " << parameter.maximum_steps_ << std::endl;
    std::cout << "thermodynamic_averaging_steps: " << parameter.thermodynamic_averaging_steps_
              << std::endl;
    std::cout << "temperature: " << parameter.temperature_ << std::endl;
    std::cout << "restart_steps: " << parameter.restart_steps_ << std::endl;
    std::cout << "restart_energy: " << parameter.restart_energy_ << std::endl;
    // std::cout << "initial_temperature: " << parameter.initial_temperature_ << std::endl;
    // std::cout << "decrement_temperature: " << parameter.decrement_temperature_ << std::endl;
  } 
  
  else if (parameter.method == "Ansys" || parameter.method == "Reformat") {
    std::cout << "initial_steps: " << parameter.initial_steps_ << std::endl;
    std::cout << "increment_steps: " << parameter.increment_steps_ << std::endl;
    std::cout << "cutoffs: ";
    std::transform(parameter.cutoffs_.begin(), parameter.cutoffs_.end(),
                   std::ostream_iterator<std::string>(std::cout, " "),
                   [](auto cutoff) { return std::to_string(cutoff); });
    std::cout << std::endl;
    std::cout << "log_type: " << parameter.log_type_ << std::endl;
    std::cout << "config_type: " << parameter.config_type_ << std::endl;
  }
}


// void Run(const Parameter &parameter) {
//   if (parameter.method == "KineticMcFirstMpi") {
//     auto kinetic_mc_first_mpi = api::BuildKineticMcFirstMpiFromParameter(parameter);
//     kinetic_mc_first_mpi.Simulate();
//   } else if (parameter.method == "KineticMcFirstOmp") {
//     auto kinetic_mc_first_omp = api::BuildKineticMcFirstOmpFromParameter(parameter);
//     kinetic_mc_first_omp.Simulate();
//   } else if (parameter.method == "KineticMcChainOmpi") {
//     auto kinetic_mc_chain_ompi = api::BuildKineticMcChainOmpiFromParameter(parameter);
//     kinetic_mc_chain_ompi.Simulate();
//   } else if (parameter.method == "Ansys") {
//     auto iterator = api::BuildIteratorFromParameter(parameter);
//     iterator.RunAnsys();
//   } else if (parameter.method == "Reformat") {
//     auto iterator = api::BuildIteratorFromParameter(parameter);
//     iterator.RunReformat();
//   } else if (parameter.method == "SimulatedAnnealing") {
//     auto simulated_annealing = api::BuildSimulatedAnnealingFromParameter(parameter);
//     simulated_annealing.Simulate();
//   } else if (parameter.method == "CanonicalMcSerial") {
//     auto canonical_mc_serial = api::BuildCanonicalMcSerialFromParameter(parameter);
//     canonical_mc_serial.Simulate();
//   } else if (parameter.method == "CanonicalMcOmp") {
//     auto canonical_mc_omp = api::BuildCanonicalMcOmpFromParameter(parameter);
//     canonical_mc_omp.Simulate();
//   } else {
//     std::cout << "No such method: " << parameter.method << std::endl;
//   }
// }

void Run(const Parameter &parameter) {
  if (parameter.method == "CanonicalMcSerial") {
    auto canonical_mc_serial = api::BuildCanonicalMcSerialFromParameter(parameter);
    canonical_mc_serial.Simulate();
  } else if (parameter.method == "KineticMcChainOmpi") {
    auto kinetic_mc_chain_ompi = api::BuildKineticMcChainOmpiFromParameter(parameter);
    // kinetic_mc_chain_ompi.Simulate();
  } else if (parameter.method == "KineticMcFirstmpi") {
    auto kinetic_mc_first_mpi = api::BuildKineticMcFirstMpiFromParameter(parameter);
    kinetic_mc_first_mpi.Simulate();
  } else if (parameter.method == "Ansys") {
    auto iterator = api::BuildIteratorFromParameter(parameter);
    iterator.RunAnsys();
  }
}

mc::CanonicalMcSerial BuildCanonicalMcSerialFromParameter(const Parameter
                                                                &parameter) {
  Config config;
  if (parameter.map_filename_.empty()) {
    // Generalized function to read configuration
    // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
    config = Config::ReadConfig(parameter.config_filename_);

    config.UpdateNeighborList(parameter.cutoffs_);

  } else {
    // Need to figure this out
    config = Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
  }

  auto supercell_config = Config::GenerateSupercell(parameter.supercell_size_,
                                                    parameter.lattice_param_,
                                                    "X",
                                                    parameter.structure_type_);
  
  supercell_config.UpdateNeighborList(parameter.cutoffs_);
  
  // Getting element set from the configuration
  auto atom_vector = config.GetAtomVector();
  std::set<Element> element_set(atom_vector.begin(),atom_vector.end());

  // Print the elements in the set
  std::cout << "element_set: ";
  for (const auto& element : element_set) {
      std::cout << element.GetElementString() << " ";
  }
  std::cout << std::endl;

  std::cout << "Finish config reading. Start CMC." << std::endl;

  
  return mc::CanonicalMcSerial{config,
                               supercell_config,
                               parameter.log_dump_steps_,
                               parameter.config_dump_steps_,
                               parameter.maximum_steps_,
                               parameter.thermodynamic_averaging_steps_,
                               parameter.restart_steps_,
                               parameter.restart_energy_,
                               parameter.temperature_,
      // parameter.initial_temperature_,
      // parameter.decrement_temperature_,
                               element_set,
                               parameter.max_cluster_size_,
                               parameter.max_bond_order_,
                               parameter.json_coefficients_filename_};
}

mc::KineticMcChainOmpi BuildKineticMcChainOmpiFromParameter(const Parameter 
                                                                  &parameter) {

  Config config;
  if (parameter.map_filename_.empty()) {
    // Generalized function to read configuration
    // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
    config = Config::ReadConfig(parameter.config_filename_);

    config.UpdateNeighborList(parameter.cutoffs_);

  } else {
    // Need to figure this out
    config = Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
  }

  auto supercell_config = Config::GenerateSupercell(parameter.supercell_size_,
                                                    parameter.lattice_param_,
                                                    "X",
                                                    parameter.structure_type_);
  
  supercell_config.UpdateNeighborList(parameter.cutoffs_);
  
  // Getting element set from the configuration
  auto atom_vector = config.GetAtomVector();
  std::set<Element> element_set(atom_vector.begin(),atom_vector.end());

  // Print the elements in the set
  std::cout << "element_set: ";
  for (const auto& element : element_set) {
      std::cout << element.GetElementString() << " ";
  }
  std::cout << std::endl;
  
  std::cout << "Finish config reading. Start KMC." << std::endl;


  return mc::KineticMcChainOmpi{config,
                                supercell_config,
                                parameter.log_dump_steps_,
                                parameter.config_dump_steps_,
                                parameter.maximum_steps_,
                                parameter.thermodynamic_averaging_steps_,
                                parameter.restart_steps_,
                                parameter.restart_energy_,
                                parameter.restart_time_,
                                parameter.temperature_,
                                element_set,
                                parameter.max_cluster_size_,
                                parameter.max_bond_order_,
                                parameter.json_coefficients_filename_,
                                parameter.time_temperature_filename_,
                                parameter.rate_corrector_,
                                parameter.vacancy_trajectory_};

}

mc::KineticMcFirstMpi BuildKineticMcFirstMpiFromParameter (const Parameter 
                                                                 &parameter) {

  Config config;
  if (parameter.map_filename_.empty()) {
    // Generalized function to read configuration
    // Supported formats are: .cfg, .POSCAR, .cfg.gz, .cfg.bz2, .POSCAR.gz, .POSCAR.bz2
    config = Config::ReadConfig(parameter.config_filename_);

    config.UpdateNeighborList(parameter.cutoffs_);

  } else {
    // Need to figure this out
    config = Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
  }

  auto supercell_config = Config::GenerateSupercell(parameter.supercell_size_,
                                                    parameter.lattice_param_,
                                                    "X",
                                                    parameter.structure_type_);
  
  supercell_config.UpdateNeighborList(parameter.cutoffs_);
  
  // Getting element set from the configuration
  auto atom_vector = config.GetAtomVector();
  std::set<Element> element_set(atom_vector.begin(),atom_vector.end());

  // Print the elements in the set
  std::cout << "element_set: ";
  for (const auto& element : element_set) {
      std::cout << element.GetElementString() << " ";
  }
  std::cout << std::endl;
  
  std::cout << "Finish config reading. Start KMC." << std::endl;

  return mc::KineticMcFirstMpi{config,
                               supercell_config,
                               parameter.log_dump_steps_,
                               parameter.config_dump_steps_,
                               parameter.maximum_steps_,
                               parameter.thermodynamic_averaging_steps_,
                               parameter.restart_steps_,
                               parameter.restart_energy_,
                               parameter.restart_time_,
                               parameter.temperature_,
                               element_set,
                               parameter.max_cluster_size_,
                               parameter.max_bond_order_,
                               parameter.json_coefficients_filename_,
                               parameter.time_temperature_filename_,
                               parameter.rate_corrector_,
                               parameter.vacancy_trajectory_};
}

ansys::Traverse BuildIteratorFromParameter(const Parameter &parameter) {

  return ansys::Traverse{parameter.initial_steps_, 
                         parameter.increment_steps_,
                         parameter.cutoffs_,
                         parameter.log_type_,
                         parameter.config_type_};
}

// mc::KineticMcFirstMpi BuildKineticMcFirstMpiFromParameter(const Parameter &parameter) {
//   std::set<Element> element_set;
//   for (const auto &element_string : parameter.element_set_) {
//     element_set.insert(Element(element_string));
//   }
//   cfg::Config config;
//   if (parameter.map_filename_.empty()) {
//     config = cfg::Config::ReadConfig(parameter.config_filename_);
//     config.ReassignLatticeVector();
//   } else {
//     config = cfg::Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
//   }
//   std::cout << "Finish config reading. Start KMC." << std::endl;
//   return mc::KineticMcFirstMpi{config,
//                                parameter.log_dump_steps_,
//                                parameter.config_dump_steps_,
//                                parameter.maximum_steps_,
//                                parameter.thermodynamic_averaging_steps_,
//                                parameter.restart_steps_,
//                                parameter.restart_energy_,
//                                parameter.restart_time_,
//                                parameter.temperature_,
//                                element_set,
//                                parameter.json_coefficients_filename_,
//                                parameter.time_temperature_filename_,
//                                parameter.rate_corrector_};
// }


// mc::KineticMcChainOmpi BuildKineticMcChainOmpiFromParameter(const Parameter &parameter) {
//   std::set<Element> element_set;
//   for (const auto &element_string : parameter.element_set_) {
//     element_set.insert(Element(element_string));
//   }
//   cfg::Config config;
//   if (parameter.map_filename_.empty()) {
//     config = cfg::Config::ReadConfig(parameter.config_filename_);
//     config.ReassignLatticeVector();
//   } else {
//     config = cfg::Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
//   }
//   std::cout << "Finish config reading. Start KMC." << std::endl;
//   return mc::KineticMcChainOmpi{config,
//                                 parameter.log_dump_steps_,
//                                 parameter.config_dump_steps_,
//                                 parameter.maximum_steps_,
//                                 parameter.thermodynamic_averaging_steps_,
//                                 parameter.restart_steps_,
//                                 parameter.restart_energy_,
//                                 parameter.restart_time_,
//                                 parameter.temperature_,
//                                 element_set,
//                                 parameter.json_coefficients_filename_,
//                                 parameter.time_temperature_filename_,
//                                 parameter.rate_corrector_};
// }
// ansys::SimulatedAnnealing BuildSimulatedAnnealingFromParameter(const Parameter &parameter) {
//   std::map<Element, size_t> solute_atom_count;
//   for (size_t i = 0; i < parameter.solute_element_set_.size(); ++i) {
//     solute_atom_count.insert(std::make_pair(Element(parameter.solute_element_set_[i]),
//                                             parameter.solute_number_set_[i]));
//   }
// 
//   return ansys::SimulatedAnnealing{{parameter.factor_, parameter.factor_, parameter.factor_},
//                                    Element(parameter.solvent_element_),
//                                    solute_atom_count,
//                                    parameter.log_dump_steps_,
//                                    parameter.config_dump_steps_,
//                                    parameter.maximum_steps_,
//                                    parameter.early_stop_steps_,
//                                    parameter.initial_temperature_,
//                                    parameter.json_coefficients_filename_};
// }
// 
// mc::CanonicalMcOmp BuildCanonicalMcOmpFromParameter(const Parameter &parameter) {
//   std::set<Element> element_set;
//   for (const auto &element_string : parameter.element_set_) {
//     element_set.insert(Element(element_string));
//   }
//   cfg::Config config;
//   if (parameter.map_filename_.empty()) {
//     config = cfg::Config::ReadConfig(parameter.config_filename_);
//     config.ReassignLatticeVector();
//   } else {
//     config = cfg::Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
//   }
//   std::cout << "Finish config reading. Start CMC." << std::endl;
//   return mc::CanonicalMcOmp{config,
//                             parameter.log_dump_steps_,
//                             parameter.config_dump_steps_,
//                             parameter.maximum_steps_,
//                             parameter.thermodynamic_averaging_steps_,
//                             parameter.restart_steps_,
//                             parameter.restart_energy_,
//                             parameter.temperature_,
//       // parameter.initial_temperature_,
//       // parameter.decrement_temperature_,
//                             element_set,
//                             parameter.json_coefficients_filename_};
// }
// 
// ansys::Traverse BuildIteratorFromParameter(const Parameter &parameter) {
//   Element vacancy("X");
//   element_set.erase(vacancy);
// 
//   return ansys::Traverse{parameter.initial_steps_, parameter.increment_steps_,
//                          Element(parameter.solvent_element_), element_set,
//                          parameter.smallest_cluster_criteria_,
//                          parameter.solvent_bond_criteria_,
//                          parameter.json_coefficients_filename_,
//                          parameter.log_type_,
//                          parameter.config_type_};
// }

} // api