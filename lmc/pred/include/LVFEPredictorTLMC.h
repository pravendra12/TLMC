#ifndef LMC_PRED_INCLUDE_LVFEPREDICTORTLMC_H_
#define LMC_PRED_INCLUDE_LVFEPREDICTORTLMC_H_

#include "Config.h"
#include "BasisSet.h"
#include "BasisSet.h"
#include "PrintUtility.h"
#include "SymmetrySpglib.h"
#include "ClusterExpansion.h"
#include "CorrelationVector.h"
#include "GetEquivalentClusters.h"
#include "ClusterExpansionParameters.h"
#include "TiledSupercell.h"
#include "LatticeSiteMapping.hpp"
#include "LatticeSiteEncodedMapping.hpp"

using namespace std;
using namespace Eigen;

class LVFEPredictorTLMC
{
public:
  LVFEPredictorTLMC(
      const ClusterExpansionParameters &ceParams,
      const TiledSupercell &tiledSupercell);

  double GetEffectiveVacancyFormationEnergy(
      const TiledSupercell &tiledSupercell,
      const LatticeSiteMapping &vacancyLatticeSite) const;

  // LatticeId to which the element will be assigned
  // Element which need to be assigned to the lattice Id
  double GetEffectiveVacancyFormationEnergy(
      const TiledSupercell &tiledSupercell,
      const LatticeSiteMapping &vacancyLatticeSite,
      const LatticeSiteMapping &targetLatticeSite,
      const Element &elementToAssign) const;

  double GetDeForVacancyMigration(
      const TiledSupercell &tiledSupercell,
      const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair) const;

private:
  VectorXd GetLocalSiteClusterVector(
      const TiledSupercell &tiledSupercell,
      const LatticeSiteMapping &vacancyLatticeSite) const;

  const size_t maxBondOrder_{};
  const size_t maxClusterSize_{};
  mutable BasisSet atomicBasis_;
  const vector<pair<vector<vector<size_t>>, LatticeClusterType>> encodedOrbitsForSite_{};
  const vector<double> lvfeECIs_{};

  /**
   * @brief Stores the symmetrically sorted neighbouring lattice IDs around each lattice site.
   *
   * This data structure is a vector of vectors, where the outer vector is indexed by the lattice site ID.
   * For each lattice site, the inner vector contains the symmetrically sorted lattice IDs that are associated with that site.
   * This mapping allows efficient lookup and processing of lattice sites and their symmetrically representations,
   * which is useful for symmetry operations and cluster vector calculations.
   */
  vector<vector<LatticeSiteEncodedMapping>> symmetricLatticeSiteEncodedMappings_{};

  void GetSymmetricallySortedLatticeIdsVectorMap(
      const TiledSupercell &tiledSupercell);

  void PrintLVFEPredictor() const;
};

#endif // LMC_PRED_INCLUDE_LVFEPREDICTORTLMC_H_
