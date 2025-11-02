#ifndef LMC_ANSYS_INCLUDE_ORDERPARAMETER_H_
#define LMC_ANSYS_INCLUDE_ORDERPARAMETER_H_

#include "Config.h"
#include <omp.h>
#include <utility>
#include <unordered_set>
#include "TiledSupercell.h"
using namespace std;

/**
 * @class SubLatticeOccupancy
 * @brief Computes the  order parameter for a given element in a configuration.
 *
 * The  order parameter is a measure of atomic ordering in binary alloys, particularly in BCC structures.
 * It quantifies the difference in site occupancy of a given element between two distinct sublattices (alpha and beta).
 */
class SubLatticeOccupancy
{
public:
  /**
   * @brief Constructs a SubLatticeOccupancy object and initializes lattice site classifications.
   *
   * Only alpha sites are initialized in the constructor as rest of the sites
   * are beta sites.
   *
   */
  SubLatticeOccupancy(const TiledSupercell &tiledSupercell);

  /**
   * @brief Compute the  order parameter (η) for each element.
   *
   * Given the fractional occupancy of each element on the alpha and beta
   * sublattices (as returned by GetAlphaBetaSiteOccupancyAll), this function
   * calculates the  order parameter for each element:
   *
   *      η = (f_alpha - f_beta) / (f_alpha + f_beta)
   *
   * where:
   *  - f_alpha: fractional occupancy of the element on alpha sites
   *  - f_beta:  fractional occupancy of the element on beta sites
   *
   * The  order parameter quantifies the degree of ordering:
   *  - η = 1   → element fully occupies the alpha sublattice
   *  - η = -1  → element fully occupies the beta sublattice
   *  - η = 0   → random distribution between sublattices
   *
   * @param elementOccupancies  Map from Element → pair(alpha occupancy, beta occupancy)
   * @return std::map<Element, double>
   *         Map from each Element to its  order parameter η.
   *
   * Reference:
   * Santodonato, L.J., Liaw, P.K., Unocic, R.R. et al.
   * Predictive multiphase evolution in Al-containing high-entropy alloys.
   * Nat Commun 9, 4520 (2018). https://doi.org/10.1038/s41467-018-06757-2
   */

  map<Element, double> ComputeOrderParameter(
      const map<Element, pair<double, double>> &elementOccupancies) const;



  /**
   * @brief Compute alpha and beta sublattice occupancies for all elements.
   *
   * This function calculates, for each chemical element present in the system,
   * the fraction of its atoms that occupy the alpha and beta sublattices in a
   * -ordered structure (e.g., BCC-based ordering). The lattice is assumed to
   * contain an equal number of alpha and beta sites.
   *
   * @param atomIndicesVector  Vector of element atomic numbers or indices for
   *                           each lattice site (length = total number of sites).
   * @param elements           Vector of all unique Element objects present in
   *                           the system (used as map keys).
   *
   * @return std::map<Element, std::pair<double, double>>
   *         A map from each Element to a pair:
   *         (alpha site occupancy fraction, beta site occupancy fraction)
   */
  map<Element, pair<double, double>> GetAlphaBetaSiteOccupancy(
      const vector<size_t> &atomIndicesVector) const;

  /**
   * @brief Retrieves the set of lattice sites classified as alpha sites.
   *
   * Alpha sites are a subset of the lattice that are identified based on second nearest
   * neighbor relationships. These sites are used to compute the  order parameter.
   *
   * @return An unordered set containing the indices of alpha lattice sites.
   */
  const unordered_set<size_t> &GetAlphaAtomIds() const;

private:
  // Store the alpha atom Ids
  const unordered_set<size_t> alphaAtomIds_{};

  /**
   * @brief Identify and return the set of alpha-sublattice atom IDs in the tiled supercell.
   *
   * This function determines which lattice sites belong to the "alpha" sublattice
   * based on second-nearest-neighbor (2NN) connectivity in the smaller reference configuration.
   * The algorithm ensures that selected alpha sites do not share any first-nearest-neighbor (1NN)
   * connections with each other, maintaining sublattice separation typical of ordered  structures.
   *
   * Algorithm overview:
   * 1. Retrieve the small (primitive) configuration from the tiled supercell.
   * 2. Build a fast lookup map for each site's first nearest neighbors.
   * 3. Iterate over all pairs of second nearest neighbors, ensuring each pair is processed only once.
   * 4. For each site in the pair, add it to the alpha sublattice set if none of its 1NNs are already alpha sites.
   * 5. Extend the resulting set of alpha lattice sites across all replicated subcells
   *    in the tiled supercell to obtain global atom IDs.
   *
   * @param tiledSupercell Reference to the tiled supercell object containing the replicated configurations.
   * @return An unordered_set<size_t> containing atom IDs corresponding to alpha-sublattice sites
   *         across the entire tiled supercell.
   *
   * @note The function assumes the small configuration’s neighbor lists are properly initialized
   *       and that the first two shells correspond to 1NN and 2NN interactions, respectively.
   */

  unordered_set<size_t> InitializeAlphaAtomIds(const TiledSupercell &tiledSupercell);
};

#endif // LMC_ANSYS_INCLUDE_ORDERING_H_