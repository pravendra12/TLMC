#ifndef LMC_PRED_INCLUDE_SYMMETRICCEPREDICTORTLMC_H_
#define LMC_PRED_INCLUDE_SYMMETRICCEPREDICTORTLMC_H_

#include "Config.h"
#include "SymmetricCE.h"
#include "PrintUtility.h"
#include "TiledSupercell.h"
#include "ClusterExpansionParameters.h"

using namespace std;

// This class need to be passes with supercellConfig object whose neighbourLists
// are updated to max cluster size else the results will be meaningless.
class SymmetricCEPredictorTLMC
{
public:
  // Constructor for TLMC
  SymmetricCEPredictorTLMC(
      const ClusterExpansionParameters &ceParams,
      const TiledSupercell &tiledSupercell,
      const Config &primitiveConfig);

  const vector<string> &GetAllowedElements() const;

  double ComputeLocalFormationEnergyOfSite(
      const TiledSupercell &tiledSupercell,
      const LatticeSiteMapping &latticeSiteMapping);

  // Computes local formation energy for a pair of lattice sites.
  // - latticeIdPair:    (site1, site2) indices in the config
  // - latticeIdPairElements: (element1, element2) assigned to those sites
  //   where element1 ↔ site1, element2 ↔ site2
  double ComputeLocalFormationEnergyForPair(
      const TiledSupercell &tiledSupercell,
      const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeIdPair,
      const pair<Element, Element> &latticeIdPairElements);

  double GetDeSwap(
      const TiledSupercell &tiledSupercell,
      const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeIdJumpPair);

  double GetDeSwap(
      const TiledSupercell &tiledSupercell,
      const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeIdPair,
      const pair<Element, Element> &latticeIdPairElements);

private:
  // Returns total formation energy
  double GetTotalFormationEnergy(
      const vector<double> &clusterVector) const;

  // Returns total energy
  double GetTotalEnergy(
      const vector<double> &clusterVector) const;

  double GetTotalEnergy(
      const vector<double> &clusterVector,
      const unordered_map<string, int> &elementCountMap) const;

  static unordered_map<string, int> GetElementCountMap(
      const Config &supercellConfig);

  void GetSymmetricallySortedLatticeEncodedMappings(
      const TiledSupercell &tiledSupercell);

  void PrintSymmetricCEPredictorInfo() const;

  const size_t nAtoms_{};
  const vector<string> allowedElements_{};
  const vector<double> clusterCutoffs_{};
  SymmetricCE symCE_;
  const vector<vector<vector<int>>> localOrbitsEncoding_;
  const unordered_map<string, double> chemicalPotentialsMap_;
  const unordered_map<string, int> elementCountMap_{};
  const vector<double> ecis_{};

  /**
   * @brief Stores the canonical sorted lattice IDs corresponding to each <LATTICE_ID, ENCODED_SMALL_CONFIG_IDX>.
   *
   * This data structure is a vector of vectors, where the outer vector is indexed by the lattice site ID.
   * For each lattice site, the inner vector contains the canonical sorted lattice IDs that are associated with that site.
   * This mapping allows efficient lookup and processing of lattice sites and their canonical representations,
   * which is useful for symmetry operations and cluster vector calculations.
   */
  vector<vector<LatticeSiteEncodedMapping>> symmetricLatticeSiteEncodedMappings_{};
};

#endif // LMC_PRED_INCLUDE_SYMMETRICCEPREDICTORTLMC_H_
