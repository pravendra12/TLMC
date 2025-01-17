/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 9/26/23 9:18 PM                                                           *
 **************************************************************************************************/

#include "McAbstract.h"
#include <utility>
#include <chrono>
#include <mpi.h>
#include "PotentialEnergyEstimator.h"

McAbstract::McAbstract(Config config,
                       Config supercell_config,
                       const unsigned long long int log_dump_steps,
                       const unsigned long long int config_dump_steps,
                       const unsigned long long int maximum_steps,
                       const unsigned long long int thermodynamic_averaging_steps,
                       const unsigned long long int restart_steps,
                       const double restart_energy,
                       const double restart_time,
                       const double temperature,
                       const std::set<Element> &element_set,
                       const size_t max_cluster_size,
                       const size_t max_bond_order,
                       const std::string &json_coefficients_filename,
                       const std::string &log_filename)
    : config_(std::move(config)),
      log_dump_steps_(log_dump_steps),
      config_dump_steps_(config_dump_steps),
      maximum_steps_(maximum_steps),
      steps_(restart_steps),
      energy_(restart_energy),
      //absolute_energy_(PotentialEnergyEstimator(
      //    json_coefficients_filename, config, supercell_config, element_set,max_cluster_size,max_bond_order).GetEnergy(config_)),
      absolute_energy_(0),
      time_(restart_time),
      temperature_(temperature),
      beta_(1.0 / constants::kBoltzmann / temperature_),
      is_restarted_(steps_ > 0),
      // thermodynamic_averaging_(thermodynamic_averaging_steps),
      generator_(static_cast<unsigned long long int>(
                     std::chrono::system_clock::now().time_since_epoch().count())),
      unit_distribution_(0.0, 1.0),
      ofs_(log_filename, is_restarted_ ? std::ofstream::app : std::ofstream::out) {
  ofs_.precision(16);

  // std::cout << "Print I am here in McAbstract" << std::endl;

  MPI_Init(nullptr, nullptr);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank_);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size_);

  // std::cout << "Print I am here .. " << std::endl;
}
McAbstract::~McAbstract() {
  MPI_Finalize();
}

