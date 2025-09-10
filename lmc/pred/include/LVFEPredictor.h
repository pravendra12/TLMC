#ifndef LMC_PRED_INCLUDE_LVFEPREDICTOR_H_
#define LMC_PRED_INCLUDE_LVFEPREDICTOR_H_

#include "Config.h"
#include "BasisSet.h"
#include "BasisSet.h"
#include "PrintUtility.h"
#include "SymmetrySpglib.h"
#include "ClusterExpansion.h"
#include "CorrelationVector.h"
#include "GetEquivalentClusters.h"
#include "ClusterExpansionParameters.h"

using namespace std;
using namespace Eigen;

class LVFEPredictor
{
public:
  LVFEPredictor(
      const ClusterExpansionParameters &ceParams,
      const Config &config);


  double GetEffectiveVacancyFormationEnergy(
      const Config &config,
      const size_t &vacancyLatticeId);

  // LatticeId to which the element will be assigned
  // Element which need to be assigned to the lattice Id
  double GetEffectiveVacancyFormationEnergy(
      const Config &config,
      const size_t &vacancyLatticeId,
      const size_t &targetLatticeId,
      const Element &elementToAssign);

  double GetDeForVacancyMigration(
    const Config &config, 
    const pair<size_t, size_t> &latticeIdJumpPair);

private:

  VectorXd GetLocalSiteClusterVector(
      const Config &config,
      const size_t &vacancyLatticeId);
  const size_t maxBondOrder_{};
  const size_t maxClusterSize_{};

  BasisSet atomicBasis_;

  const vector<pair<vector<vector<size_t>>, LatticeClusterType>> encodedOrbitsForSite_{};
  const vector<double> lvfeECIs_{};
};

#endif // LMC_PRED_INCLUDE_LVFEPREDICTOR_H_
