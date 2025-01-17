#include "KineticMcFirstMpi.h"

// First Order Residence Time Algorithm
//    k -> i
//    |
// Current Position

namespace mc {

KineticMcFirstMpi::KineticMcFirstMpi(Config config,
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
                                     const std::string &time_temperature_filename,
                                     const bool is_rate_corrector,
                                     const Eigen::RowVector3d &vacancy_trajectory)
    : KineticMcFirstAbstract(std::move(config),
                             supercell_config,
                             log_dump_steps,
                             config_dump_steps,
                             maximum_steps,
                             thermodynamic_averaging_steps,
                             restart_steps,
                             restart_energy,
                             restart_time,
                             temperature,
                             element_set,
                             max_cluster_size,
                             max_bond_order,
                             json_coefficients_filename,
                             time_temperature_filename,
                             is_rate_corrector,
                             vacancy_trajectory) {
  if (world_size_ != kEventListSize) {
    std::cout << "Must use " << kEventListSize << " processes. Terminating...\n" << std::endl;
    MPI_Finalize();
    exit(0);
  }
  if (world_rank_ == 0) {
    std::cout << "Using " << world_size_ << " processes." << std::endl;
  }
}
KineticMcFirstMpi::~KineticMcFirstMpi() = default;

void KineticMcFirstMpi::BuildEventList() {
  total_rate_k_ = 0;
  // Neighbours of Vacancy 
  const auto neighbor_vacancy_id = 
      config_.GetNeighborLatticeIdVectorOfLattice(vacancy_lattice_id_, 1)[static_cast<size_t>(world_rank_)];

  JumpEvent event_k_i_({vacancy_lattice_id_, neighbor_vacancy_id},
      vacancy_migration_predictor_.GetBarrierAndDiffFromLatticeIdPair(
          config_, {vacancy_lattice_id_, neighbor_vacancy_id}),
      beta_);  

  const double rate_k_i = event_k_i_.GetForwardRate();


  MPI_Allreduce(&rate_k_i, &total_rate_k_, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  event_k_i_.CalculateProbability(total_rate_k_);

  MPI_Allgather(&event_k_i_,
                sizeof(JumpEvent),
                MPI_BYTE,
                event_k_i_list_.data(),
                sizeof(JumpEvent),
                MPI_BYTE,
                MPI_COMM_WORLD);
  double cumulative_probability = 0.0;
  for (auto& event_it : event_k_i_list_) {
    cumulative_probability += event_it.GetProbability();
    event_it.SetCumulativeProbability(cumulative_probability);
  }
}

double KineticMcFirstMpi::CalculateTime() {
  static std::uniform_real_distribution<double> uniform_distribution(0.0, 1.0);
  return -std::log(uniform_distribution(generator_)) / total_rate_k_ / constants::kPrefactor;
}
} // mc