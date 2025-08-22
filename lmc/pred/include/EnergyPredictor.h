#ifndef LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#define LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_

#include "Config.h"
#include "ClusterExpansion.h"
#include "SymmetrySpglib.h"
#include "JsonUtility.h"
#include <string>
#include "ClusterExpansionParameters.h"
#include "CorrelationVector.h"
#include "BasisSet.h"

using namespace std;

class EnergyPredictor
{
public:
  EnergyPredictor(
    const ClusterExpansionParameters &ceParams,
      const Config &config);

  // Returns the energy of a site
  double ComputeEnergyOfSite(
      const Config &config,
      const size_t &latticeId);

  // Returns total energy of config using Symmetric site centered cluster expansion
  // E = sum_i * Ei
  double ComputeEnergyOfConfig(const Config &config);

  // Returns the dE due to vacancy migration
  double GetDeMigration(
    const Config &config, 
    const pair<size_t, size_t> &latticeIdJumpPair);

  // Returns the dE due to atom swap 
  double GetDeSwap(
    Config &config, 
    const pair<size_t, size_t> &latticeIdJumpPair);

  // Returns the energy change due to atom swap
  double ComputeDeltaEnergy(
      const Config &config,
      const pair<size_t, size_t> &latticeIdPair);

private:
  /*! @brief Maximum Cluster Size
   */
  const size_t maxClusterSize_{};

  /*! @brief Maximum Bond Order
   */
  const size_t maxBondOrder_{};

  /*! @brief Effective cluster interactions
   */
  const VectorXd ecis_{};

  /*! @brief Basis Set
  */
  BasisSet atomicBasis_;

  /*! @brief Element set
   */
  const set<Element> elementSet_{};

  /*! @brief Equivalent Clusters Encoding
   */
  const vector<pair<vector<vector<size_t>>,
                    LatticeClusterType>>
      equivalentClustersEncoding_;
};

#endif // LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_

