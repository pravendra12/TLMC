#ifndef LMC_PRED_INCLUDE_KRAPREDICTOR_H_
#define LMC_PRED_INCLUDE_KRAPREDICTOR_H_

#include <Eigen/Dense>
#include <string>
#include "Config.h"
#include "Symmetry.h"
#include "GetOrbits.h"
#include "JsonUtility.h"
#include "LocalEnvironmentEncoder.h"

using namespace std;
using namespace Eigen;

class KRAPredictor
{
public:
  KRAPredictor(const string &predictorFilename,
               const Config &config,
               const set<Element> &elementSet,
               const size_t maxBondOrder,
               const size_t maxClusterSize);

  // Returns Kinetically resolved activation barrier
  [[nodiscard]] double GetKRA(
      const Config &config,
      const pair<size_t, size_t> &latticeIdJumpPair) const;

private:
  // Need to think about to make it more generalized so that based on the initialization
  // of the simulation it can be declared

  // Fitting coefficients
  const VectorXd betaKRA_W_;
  const double interceptKRA_W_;

  const VectorXd betaKRA_Ta_;
  const double interceptKRA_Ta_;

  // Element set (without X)
  const set<Element> elementSet_;

  // CE details
  const size_t maxBondOrder_;
  const size_t maxClusterSize_;

  // Atomic Basis
  const string basisType_ = "Occupation";

  // Equivalent sites encoding under 3 Bar symmetry
  vector<vector<size_t>> equivalentSites3Bar_;
};

#endif // LMC_PRED_INCLUDE_KRAPREDICTOR_H_