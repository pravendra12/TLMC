#ifndef LMC_CE_INCLUDE_LATTICECLUSTERMMM_HPP_
#define LMC_CE_INCLUDE_LATTICECLUSTERMMM_HPP_

#include <iostream>
#include "Config.h"
#include "LatticeCluster.hpp"
#include "LatticeClusterType.hpp"




/*! \brief Comparator for sorting lattice site pairs based on their positions.
 *  \details 
 *  - Adjusts positions by subtracting 0.5 from each component.
 *  - Compares the squared norm of adjusted positions as the primary criterion.
 *  - Resolves ties using:
 *    1. Absolute x-symmetry.
 *    2. Individual x, y, and z components in sequence.
 *  \param lhs : First pair of lattice ID and position.
 *  \param rhs : Second pair of lattice ID and position.
 *  \return True if lhs should precede rhs in sorting; otherwise, false.
 */

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

// singlets
inline bool IsClusterSmallerSymmetricallyMMMSinglets(const std::pair<size_t, Eigen::RowVector3d>& lhs,
                                             const std::pair<size_t, Eigen::RowVector3d>& rhs) {
   
  if (GroupCompareMMM(lhs, rhs)) {
    return true;
  }
  if (GroupCompareMMM(rhs, lhs)) {
    return false;
  }

  return false;

}
// for pairs and triplets

inline bool IsClusterSmallerSymmetricallyMMMSinglets(const std::vector<std::pair<size_t, Eigen::RowVector3d>>& lhs,
                                                     const std::vector<std::pair<size_t, Eigen::RowVector3d>>& rhs) {
  for (size_t i = 0; i<lhs.size(); i++){ 

    const auto& lhs_lattice = lhs[i];
    const auto& rhs_lattice = rhs[i];

    if (GroupCompareMMM(lhs_lattice, rhs_lattice)) {
      return true;
    }
    if (GroupCompareMMM(rhs_lattice, lhs_lattice)) {
      return false;
    }
  }

  return false;

}


/*! \brief Class for defining an extended LatticeClusterMMM
*/

class LatticeClusterMMM {
 public:
  /*! \brief Default constructor for LatticeClusterMMM.
   */
  LatticeClusterMMM() = default;

  /*! \brief Constructor for setting up the extended cluster of lattice sites.
   *  \param lattice_cluster_type : The cluster type.
   *  \param lattice_id_vector    : The lattice id vector of the cluster, which 
   *                                contains pair of lattice id and position.
   *  \param property_vector      : An additional property vector for the cluster.
   */
  LatticeClusterMMM(std::vector<std::pair<size_t, Eigen::RowVector3d>> lattice_id_vector)
      : lattice_id_vector_(std::move(lattice_id_vector)) {
       
       Sort();
       symmetry_label_ = FindSymmetryLabel();

      }

  /*! \brief Default destructor for LatticeClusterMMM.
   */
  virtual ~LatticeClusterMMM() = default;

  void checkFunc(){
    std::cout << "Yes it is working" << std::endl;
  }

    void Sort() {
        std::sort(lattice_id_vector_.begin(), lattice_id_vector_.end(), PositionCompareMMM);
    }

    bool FindSymmetryLabel() {
    for (size_t i = 0; i < lattice_id_vector_.size(); ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (!GroupCompareMMM(lattice_id_vector_[i], lattice_id_vector_[j]) &&
                !GroupCompareMMM(lattice_id_vector_[j], lattice_id_vector_[i])) {
                return true;
            }
        }
    }
    return false;
}


  private:

    /// Stores lattice id and relative position as a pair.
    std::vector<std::pair<size_t, Eigen::RowVector3d>> lattice_id_vector_;
    bool symmetry_label_;

  

};
#endif //LMC_CE_INCLUDE_LATTICECLUSTER_HPP_
