#ifndef LMC_PRED_INCLUDE_SYMMETRICCEPREDICTOR_H_
#define LMC_PRED_INCLUDE_SYMMETRICCEPREDICTOR_H_

#include "Config.h"
#include "SymmetricCE.h"
#include "PrintUtility.h"
#include "ClusterExpansionParameters.h"

using namespace std;

// This class need to be passes with supercellConfig object whose neighbourLists
// are updated to max cluster size else the results will be meaningless.
class SymmetricCEPredictor
{
public:
  SymmetricCEPredictor(
      const ClusterExpansionParameters &ceParams,
      const Config &supercellConfig,
      const Config &primitiveConfig);

  const vector<string> &GetAllowedElements() const;

  double ComputeLocalFormationEnergyOfSite(
      const Config &config,
      const size_t &latticeId);

  // Returns total energy of config using Symmetric site centered cluster expansion
  // E = sum_i * Ei
  double ComputeEnergyOfConfig(const Config &config);

  // Returns the dE due to atom swap
  double GetDeSwap(
      Config &config,
      const pair<size_t, size_t> &latticeIdJumpPair);

  // Computes local formation energy for a pair of lattice sites.
  // - latticeIdPair:    (site1, site2) indices in the config
  // - latticeIdPairElements: (element1, element2) assigned to those sites
  //   where element1 ↔ site1, element2 ↔ site2
  double ComputeLocalFormationEnergyForPair(
      const Config &config,
      const pair<size_t, size_t> &latticeIdPair,
      const pair<Element, Element> &latticeIdPairElements);

  // Need to rename it
  double GetDeSwapConst(
      const Config &config,
      const pair<size_t, size_t> &latticeIdJumpPair);

  double GetDeSwapConst(
      const Config &config,
      const pair<size_t, size_t> &latticeIdPair,
      const pair<Element, Element> &latticeIdPairElements);

private:
  // Returns total formation energy
  double GetTotalFormationEnergy(
      const vector<double> &clusterVector);

  // Returns total energy
  double GetTotalEnergy(
      const vector<double> &clusterVector);

  double GetTotalEnergy(
      const vector<double> &clusterVector,
      const unordered_map<string, int> &elementCountMap);

  static unordered_map<string, int> GetElementCountMap(
      const Config &supercellConfig);

  static vector<vector<size_t>> GetSymmetricallySortedLatticeIdsVectorMap(
      const Config &supercellConfig);

  const size_t nAtoms_{};
  const vector<string> allowedElements_{};
  const vector<double> clusterCutoffs_{};
  SymmetricCE symCE_;
  const vector<vector<vector<int>>> localOrbitsEncoding_;
  const unordered_map<string, double> chemicalPotentialsMap_;
  const unordered_map<string, int> elementCountMap_{};
  const vector<double> ecis_{};

  /**
   * @brief Stores the canonical sorted lattice IDs corresponding to each lattice site.
   *
   * This data structure is a vector of vectors, where the outer vector is indexed by the lattice site ID.
   * For each lattice site, the inner vector contains the canonical sorted lattice IDs that are associated with that site.
   * This mapping allows efficient lookup and processing of lattice sites and their canonical representations,
   * which is useful for symmetry operations and cluster vector calculations.
   */
  const vector<vector<size_t>> symmetricallySortedLatticeIdsVectorMap_{};
};

#endif // LMC_PRED_INCLUDE_SYMMETRICCEPREDICTOR_H_
