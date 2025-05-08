#ifndef LMC_CE_INCLUDE_SYMMETRY_H_
#define LMC_CE_INCLUDE_SYMMETRY_H_

#include "Config.h"
#include "Eigen/Dense"
#include "Constants.hpp"
#include "EncodingUtility.h"
#include "PrintUtility.h"

using namespace std;
using namespace Eigen;


struct Vector3dHash
{
  std::size_t operator()(const Vector3d &v) const
  {
    return std::hash<double>()(v.x()) ^ std::hash<double>()(v.y()) ^ std::hash<double>()(v.z());
  }
};

/**
 * @brief Computes and returns a sorted list of lattice site IDs for a given pair of lattice
 *        sites under 3-bar symmetry. The sorting is based on the distance from the center of
 *        the pair and their relative positions after applying a symmetry transformation.
 *
 * @param config The configuration object containing lattice information and utility functions.
 * @param latticeIdJumpPair A pair of lattice site IDs representing the lattice pair.
 * @param maxBondOrder The maximum bond order to consider for neighboring lattice sites.
 * @return A vector of lattice site IDs sorted based on their transformed positions and distances.
 *
 * The function performs the following steps:
 * 1. Retrieves the neighboring lattice site IDs for the given pair up to the specified bond
 *    order.
 * 2. Computes the center of the lattice pair and determines the rotation axis based on the
 *    pair's direction and the closest [111] crystallographic direction.
 * 3. Constructs a rotation matrix to align the rotation axis with the z-axis.
 * 4. Transforms the positions of the neighboring lattice sites using the rotation matrix and
 *    computes a sorting key based on their distances from the center and their relative
 *    positions.
 * 5. Sorts the lattice site IDs based on the computed keys.
 * 6. Returns the sorted list of lattice site IDs.
 *
 * The sorting ensures that lattice sites are ordered consistently under the 3-bar symmetry
 * transformation, which is useful for symmetry-aware computations in lattice-based simulations.
 */
vector<size_t> GetSortedLatticeStatesForPairUnder3BarSymmetry(
  const Config &config, 
  const pair<size_t, size_t> &latticeIdJumpPair, 
  const size_t &maxBondOrder);


  /**
 * @brief Computes the equivalent lattice sites under 3-bar symmetry for a given lattice pair.
 *
 * This function determines the equivalent lattice sites under 3-bar symmetry
 * for a given lattice pair in the configuration. It uses the symmetry properties
 * to group lattice sites into equivalence classes based on their positions and
 * symmetry transformations.
 *
 * @param config The configuration object containing lattice information.
 * @param maxBondOrder The maximum bond order to consider for symmetry calculations.
 * @param latticeIdPair A pair of lattice IDs representing the central and neighboring lattice sites.
 * @return A vector of vectors, where each inner vector contains the indices of equivalent lattice sites.
 *
 * The function performs the following steps:
 * 1. Retrieves the sorted lattice vector state for the given pair under 3-bar symmetry.
 * 2. Computes the rotation axis and center of the lattice pair.
 * 3. Maps Cartesian positions of lattice sites to their indices.
 * 4. Iterates through the lattice sites, identifying equivalent sites using symmetry transformations.
 * 5. Groups equivalent sites into equivalence classes and sorts them.
 *
 * Note:
 * - The function assumes that the configuration object provides methods to retrieve
 *   Cartesian positions of lattice sites and other necessary information.
 * - The symmetry transformations are applied using helper functions such as
 *   `GetEquivalentPositionUnder3BarSymmetry`.
 */
vector<vector<size_t>> GetEquivalentSitesUnder3BarSymmetry(
    const Config &config,
    size_t maxBondOrder,
    const pair<size_t, size_t> &latticeIdPair);

#endif // LMC_CE_INCLUDE_SYMMETRY_H_
