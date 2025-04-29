#include "CanonicalMcOmp.h"
#include <omp.h>
namespace mc
{
  CanonicalMcOmp::CanonicalMcOmp(Config config,
                                 Config supercell_config,
                                 unsigned long long int log_dump_steps,
                                 unsigned long long int config_dump_steps,
                                 unsigned long long int maximum_steps,
                                 unsigned long long int thermodynamic_averaging_steps,
                                 unsigned long long int restart_steps,
                                 double restart_energy,
                                 double temperature,
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
                            element_set,
                            max_cluster_size,
                            max_bond_order,
                            json_coefficients_filename)
  {
    if (world_size_ != 1)
    {
      std::cout << "Must use 1 process. Terminating...\n"
                << std::endl;
      MPI_Finalize();
      exit(0);
    }
#pragma omp parallel default(none) shared(std::cout)
    {
#pragma omp master
      {
        std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
        num_threads_ = static_cast<size_t>(omp_get_num_threads());
      }
    }
    event_vector_.reserve(static_cast<size_t>(omp_get_num_threads()));
  }
  void CanonicalMcOmp::BuildEventVector()
  {
    event_vector_.clear();
    unavailable_position_.clear();

    std::pair<size_t, size_t> lattice_id_jump_pair;
    for (size_t i = 0; i < num_threads_; ++i)
    {
      size_t ct = 0;
      do
      {
        lattice_id_jump_pair = GenerateLatticeIdJumpPair();
        ct++;
        if (ct == 50)
        {
          break;
        }
      } while (unavailable_position_.find(lattice_id_jump_pair.first) != unavailable_position_.end() || unavailable_position_.find(lattice_id_jump_pair.second) != unavailable_position_.end());
      if (ct == 50)
      {
        break;
      }
      for (auto selected_lattice_index : {lattice_id_jump_pair.first, lattice_id_jump_pair.second})
      {
        unavailable_position_.emplace(selected_lattice_index);

        auto firstNN = config_.GetNeighborLatticeIdVectorOfLattice(selected_lattice_index, 1);
        std::copy(firstNN.begin(), firstNN.end(),
                  std::inserter(unavailable_position_,
                                unavailable_position_.end()));

        auto secondNN = config_.GetNeighborLatticeIdVectorOfLattice(selected_lattice_index, 2);
        std::copy(secondNN.begin(), secondNN.end(),
                  std::inserter(unavailable_position_,
                                unavailable_position_.end()));

        auto thirdNN = config_.GetNeighborLatticeIdVectorOfLattice(selected_lattice_index, 3);
        std::copy(thirdNN.begin(), thirdNN.end(),
                  std::inserter(unavailable_position_,
                                unavailable_position_.end()));
      }
      event_vector_.emplace_back(lattice_id_jump_pair, 0);
    }
#pragma omp parallel for default(none)
    for (auto &event : event_vector_)
    {
      event.second = energy_change_predictor_.GetDeThreadSafe(config_, event.first);
    }
  }
  void CanonicalMcOmp::Simulate()
  {
    // while (steps_ <= maximum_steps_ * static_cast<unsigned long long int>(initial_temperature_ / decrement_temperature_ + 1)) {
    while (steps_ <= maximum_steps_)
    {
      BuildEventVector();
      for (auto [lattice_id_jump_pair, dE] : event_vector_)
      {
        thermodynamic_averaging_.AddEnergy(energy_);
        Dump();
        // UpdateTemperature();
        SelectEvent(lattice_id_jump_pair, dE);
        ++steps_;
      }
    }
  }
} // mc