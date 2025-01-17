#include "SymmetryCustom.h"


std::unordered_set<size_t> 
GetNeighboringLatticeIdsOfPair(const Config& config, 
                               const std::pair<size_t, size_t>& lattice_id_pair, 
                               const size_t& max_bond_order) {

  // unordered set to store the neighboring_lattice_ids 
  std::unordered_set<size_t> neighboring_lattice_ids;

  // Add neighbors for lattice_id1 up to the 3rd nearest neighbor, excluding lattice id Jump Pair
  for (int bond_order = 1; bond_order <= max_bond_order; ++bond_order) {
    
    auto neighbors_vector = config.GetNeighborLatticeIdVectorOfLattice(lattice_id_pair.first, bond_order);
    
    for (const auto& neighbor_id : neighbors_vector) {
      if (neighbor_id != lattice_id_pair.first && neighbor_id != lattice_id_pair.second) {
        
        neighboring_lattice_ids.insert(neighbor_id);
      }
    }
  }
  // Add neighbors for lattice_id2 up to the 3rd nearest neighbor, excluding lattice_id1 and lattice_id2
  for (int bond_order = 1; bond_order <= max_bond_order; ++bond_order) {
    
    auto neighbors_vector = config.GetNeighborLatticeIdVectorOfLattice(lattice_id_pair.second, bond_order);
    
    for (const auto& neighbor_id : neighbors_vector) {
      if (neighbor_id != lattice_id_pair.first && neighbor_id != lattice_id_pair.second) {
        
        neighboring_lattice_ids.insert(neighbor_id);
      }
    }
  }
  
  return neighboring_lattice_ids;
}