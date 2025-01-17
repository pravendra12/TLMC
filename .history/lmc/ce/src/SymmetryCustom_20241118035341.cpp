#include "SymmetryCustom.h"

std::unordered_set<size_t> 
GetNeighboringLatticeIdsOfPair(const  Config& config, 
                               const std::pair<size_t, size_t>& lattice_id_pair
                               const size_t& max_bond_order) {
  // unordered set to store the neigbouring lattice ids
  std::unordered_set<size_t> neighboring_lattice_ids;

  // Add neighbors for lattice_id1 up to the 3rd nearest neighbor, excluding lattice id Jump Pair
  for (int nn_order = 1; nn_order <= 3; ++nn_order) {
      auto neighbors = cfg.GetNeighborLatticeIdVectorOfLattice(lattice_id1, nn_order);
      for (const auto& neighbor : neighbors) {
          if (neighbor != lattice_id1 && neighbor != lattice_id2) {
              unique_neighbors.insert(neighbor);
          }
      }
  }
  // Add neighbors for lattice_id2 up to the 3rd nearest neighbor, excluding lattice_id1 and lattice_id2
  for (int nn_order = 1; nn_order <= 3; ++nn_order) {
      auto neighbors = cfg.GetNeighborLatticeIdVectorOfLattice(lattice_id2, nn_order);
      for (const auto& neighbor : neighbors) {
          if (neighbor != lattice_id1 && neighbor != lattice_id2) {
              unique_neighbors.insert(neighbor);
          }
      }
  }
  return unique_neighbors;
}