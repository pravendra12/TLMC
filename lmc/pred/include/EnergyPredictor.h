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
#include "SymmetricCE.h"

using namespace std;

class EnergyPredictor
{
public:
  EnergyPredictor(
      const string &predictorFilename,
      const Config &supercellConfig,
      const Config &primitiveConfig,
      const vector<string> &allowedElements,
      const vector<double> &clusterCutoffs);

  double ComputeLocalFormationEnergyOfSite(
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
  // Returns total formation energy
  double GetTotalFormationEnergy(
      const vector<double> &clusterVector);

  // Returns total energy
  double GetTotalEnergy(
      const vector<double> &clusterVector);

  /*! @brief Effective cluster interactions
   */
  const vector<double> ecis_{};

  SymmetricCE symCE_;

  const vector<vector<vector<int>>> localOrbitsEncoding_;
  const size_t nAtoms_;

  unordered_map<string, int> elementCountMap_;

  // Chemical potentials (example: Mo, Ta)
  const unordered_map<string, double> chemicalPotentials_{
      {"Mo", -10.93308598}, // Mo (Z=42)
      {"Ta", -11.81233671}  // Ta (Z=73)
  };
};

#endif // LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
