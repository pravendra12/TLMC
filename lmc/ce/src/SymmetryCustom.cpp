#include "SymmetryCustom.h"
#include <algorithm>







void RotateLatticeVector(std::unordered_map<size_t, Eigen::RowVector3d>& lattice_id_hashmap, const Eigen::Matrix3d& rotation_matrix) {

  const Eigen::RowVector3d move_distance_after_rotation = Eigen::RowVector3d(0.5, 0.5, 0.5) - (Eigen::RowVector3d(0.5, 0.5, 0.5) * rotation_matrix);

  for (auto &lattice: lattice_id_hashmap) {

    auto relative_position = lattice.second;
    // rotate
    relative_position = relative_position * rotation_matrix;
    // move to new center
    relative_position += move_distance_after_rotation;
    relative_position -= relative_position.unaryExpr([](double x) { return std::floor(x); });

    lattice.second = relative_position;
  }
}
// 
//  
//  
//  
//  

// position compare mmm
inline bool PositionCompareMMM(const std::pair<size_t, Eigen::RowVector3d>& lhs,
                               const std::pair<size_t, Eigen::RowVector3d>& rhs) {
    const auto& relative_position_lhs = lhs.second;
    const auto& relative_position_rhs = rhs.second;

    // Adjusting the positions (subtracting 0.5 from each component of the vector)
    Eigen::RowVector3d adjusted_lhs = relative_position_lhs - Eigen::RowVector3d::Constant(0.5);
    Eigen::RowVector3d adjusted_rhs = relative_position_rhs - Eigen::RowVector3d::Constant(0.5);

    // Compute the difference in the adjusted positions
    const double diff_norm = adjusted_lhs.squaredNorm() - adjusted_rhs.squaredNorm();
    if (diff_norm < -constants::kEpsilon) {
        return true;
    }
    if (diff_norm > constants::kEpsilon) {
        return false;
    }

    // Compare the absolute x-symmetry
    const double diff_x_sym = std::abs(adjusted_lhs[0]) - std::abs(adjusted_rhs[0]);
    if (diff_x_sym < -constants::kEpsilon) {
        return true;
    }
    if (diff_x_sym > constants::kEpsilon) {
        return false;
    }

    // Compare individual components (x, y, z)
    const double diff_x = adjusted_lhs[0] - adjusted_rhs[0];
    if (diff_x < -constants::kEpsilon) {
        return true;
    }
    if (diff_x > constants::kEpsilon) {
        return false;
    }

    const double diff_y = adjusted_lhs[1] - adjusted_rhs[1];
    if (diff_y < -constants::kEpsilon) {
        return true;
    }
    if (diff_y > constants::kEpsilon) {
        return false;
    }

    const double diff_z = adjusted_lhs[2] - adjusted_rhs[2];
    if (diff_z < -constants::kEpsilon) {
        return true;
    }
    if (diff_z > constants::kEpsilon) {
        return false;
    }

    return false;
}

inline bool PositionCompareMM2(const std::pair<size_t, Eigen::RowVector3d>& lhs,
                               const std::pair<size_t, Eigen::RowVector3d>& rhs) {
    const auto& relative_position_lhs = lhs.second;
    const auto& relative_position_rhs = rhs.second;

    // Adjusting the positions (subtracting 0.5 from each component of the vector)
    Eigen::RowVector3d adjusted_lhs = relative_position_lhs - Eigen::RowVector3d::Constant(0.5);
    Eigen::RowVector3d adjusted_rhs = relative_position_rhs - Eigen::RowVector3d::Constant(0.5);

    // Compute the difference in the adjusted positions
    const double diff_norm = adjusted_lhs.squaredNorm() - adjusted_rhs.squaredNorm();
    if (diff_norm < -constants::kEpsilon) {
        return true;
    }
    if (diff_norm > constants::kEpsilon) {
        return false;
    }

    // Compare individual components (x, y, z)
    const double diff_x = relative_position_lhs[0] - relative_position_rhs[0];
    if (diff_x < -constants::kEpsilon) {
        return true;
    }
    if (diff_x > constants::kEpsilon) {
        return false;
    }

    const double diff_y = relative_position_lhs[1] - relative_position_rhs[1];
    if (diff_y < -constants::kEpsilon) {
        return true;
    }
    if (diff_y > constants::kEpsilon) {
        return false;
    }

    const double diff_z = relative_position_lhs[2] - relative_position_rhs[2];
    if (diff_z < -constants::kEpsilon) {
        return true;
    }
    if (diff_z > constants::kEpsilon) {
        return false;
    }

    return false;
}


std::vector<size_t> GetSymmetricallySortedLatticeVectorMMM(
    const Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair, 
    const size_t max_bond_order) {

  // Get neighboring lattice IDs
  auto neighboring_lattice_ids = config.GetNeighboringLatticeIdSetOfPair(lattice_id_jump_pair, max_bond_order);

  size_t num_sites = neighboring_lattice_ids.size();
  Eigen::RowVector3d move_distance = Eigen::RowVector3d(0.5, 0.5, 0.5) - config.GetLatticePairCenter(lattice_id_jump_pair);

  std::unordered_map<size_t, Eigen::RowVector3d> lattice_id_hashmap;
  lattice_id_hashmap.reserve(num_sites);

  // Move lattice IDs to center
  for (const auto id : neighboring_lattice_ids) {
    Eigen::RowVector3d relative_position = config.GetRelativePositionOfLattice(id).transpose();
    relative_position += move_distance;
    relative_position -= relative_position.unaryExpr([](double x) { return std::floor(x); });
    lattice_id_hashmap.emplace(id, relative_position);
  }

  // Rotate lattice vectors
  RotateLatticeVector(lattice_id_hashmap, config.GetLatticePairRotationMatrix(lattice_id_jump_pair));

  // Convert unordered_map to vector for sorting
  std::vector<std::pair<size_t, Eigen::RowVector3d>> lattice_id_vector(lattice_id_hashmap.begin(), lattice_id_hashmap.end());

  // Sort the lattice vector based on PositionCompareMMM
  std::sort(lattice_id_vector.begin(), lattice_id_vector.end(), PositionCompareMMM);
  
  std::cout << "MMM Symmetrically Sorted Positions : " << std::endl;
  for (const auto& pair : lattice_id_vector){
    std::cout << pair.first << " : " << pair.second << std::endl;
  }

  // Extract and return only the lattice IDs
  std::vector<size_t> sorted_lattice_ids;
  sorted_lattice_ids.reserve(lattice_id_vector.size());
  for (const auto& pair : lattice_id_vector) {
    sorted_lattice_ids.push_back(pair.first);
  }

  return sorted_lattice_ids;
}


std::vector<size_t> GetSymmetricallySortedLatticeVectorMM2(
    const Config &config, const std::pair<size_t,size_t> &lattice_id_jump_pair,
    const size_t max_bond_order){
  // Get neighboring lattice IDs
  auto neighboring_lattice_ids = config.GetNeighboringLatticeIdSetOfPair(lattice_id_jump_pair, max_bond_order);

  size_t num_sites = neighboring_lattice_ids.size();
  Eigen::RowVector3d move_distance = Eigen::RowVector3d(0.5, 0.5, 0.5) - config.GetLatticePairCenter(lattice_id_jump_pair);

  std::unordered_map<size_t, Eigen::RowVector3d> lattice_id_hashmap;
  lattice_id_hashmap.reserve(num_sites);

  // Move lattice IDs to center
  for (const auto id : neighboring_lattice_ids) {
    Eigen::RowVector3d relative_position = config.GetRelativePositionOfLattice(id).transpose();
    relative_position += move_distance;
    relative_position -= relative_position.unaryExpr([](double x) { return std::floor(x); });
    lattice_id_hashmap.emplace(id, relative_position);
  }

  // Rotate lattice vectors
  RotateLatticeVector(lattice_id_hashmap, config.GetLatticePairRotationMatrix(lattice_id_jump_pair));

  // Convert unordered_map to vector for sorting
  std::vector<std::pair<size_t, Eigen::RowVector3d>> lattice_id_vector(lattice_id_hashmap.begin(), lattice_id_hashmap.end());

  // Sort the lattice vector based on PositionCompareMMM
  std::sort(lattice_id_vector.begin(), lattice_id_vector.end(), PositionCompareMM2);
  

  std::cout << "MM2 Symmetrically Sorted Positions : " << std::endl;
  for (const auto& pair : lattice_id_vector){
    std::cout << pair.first << " : " << pair.second << std::endl;
  }

  // Extract and return only the lattice IDs
  std::vector<size_t> sorted_lattice_ids;
  sorted_lattice_ids.reserve(lattice_id_vector.size());
  for (const auto& pair : lattice_id_vector) {
    sorted_lattice_ids.push_back(pair.first);
  }

  return sorted_lattice_ids;
}




inline bool GroupCompareMMM(const std::pair<size_t, Eigen::RowVector3d>& lhs,
                            const std::pair<size_t, Eigen::RowVector3d>& rhs) {
    const auto& relative_position_lhs = lhs.second;
    const auto& relative_position_rhs = rhs.second;
     
    // Adjusting the positions (subtracting 0.5 from each component of the vector)
    Eigen::RowVector3d adjusted_lhs = relative_position_lhs - Eigen::RowVector3d::Constant(0.5);
    Eigen::RowVector3d adjusted_rhs = relative_position_rhs - Eigen::RowVector3d::Constant(0.5);

    // Compute the difference in the adjusted positions
    const double diff_norm = adjusted_lhs.squaredNorm() - adjusted_rhs.squaredNorm();
    if (diff_norm < -constants::kEpsilon) {
        return true;
    }
    if (diff_norm > constants::kEpsilon) {
        return false;
    }

    // Compare the absolute x-symmetry
    const double diff_x_sym = std::abs(adjusted_lhs[0]) - std::abs(adjusted_rhs[0]);
    return diff_x_sym < -constants::kEpsilon;
}


inline bool GroupCompareMM2(const std::pair<size_t, Eigen::RowVector3d>& lhs,
                            const std::pair<size_t, Eigen::RowVector3d>& rhs) {
    const auto& relative_position_lhs = lhs.second;
    const auto& relative_position_rhs = rhs.second;
     
    // Adjusting the positions (subtracting 0.5 from each component of the vector)
    Eigen::RowVector3d adjusted_lhs = relative_position_lhs - Eigen::RowVector3d::Constant(0.5);
    Eigen::RowVector3d adjusted_rhs = relative_position_rhs - Eigen::RowVector3d::Constant(0.5);

    // Compute the difference in the adjusted positions
    const double diff_norm = adjusted_lhs.squaredNorm() - adjusted_rhs.squaredNorm();
    if (diff_norm < -constants::kEpsilon) {
        return true;
    }
    if (diff_norm > constants::kEpsilon) {
        return false;
    }

    // Compare the absolute x-symmetry
    const double diff_x = relative_position_lhs[0] - relative_position_rhs[0];
    return diff_x < -constants::kEpsilon;
}




    