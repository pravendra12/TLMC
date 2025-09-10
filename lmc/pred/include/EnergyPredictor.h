#ifndef LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_
#define LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_

#include "Config.h"
#include "SymmetricCE.h"
#include "PrintUtility.h"
#include "ClusterExpansionParameters.h"

using namespace std;

class EnergyPredictor
{
public:
    EnergyPredictor(
        const ClusterExpansionParameters &ceParams,
        const Config &supercellConfig,
        const Config &primitiveConfig);

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


  private:
    // Returns total formation energy
    double GetTotalFormationEnergy(
        const vector<double> &clusterVector);

    // Returns total energy
    double GetTotalEnergy(
        const vector<double> &clusterVector);

    static unordered_map<string, int> GetElementCountMap(
        const Config &supercellConfig);
    
    const size_t nAtoms_{};
    const vector<string> allowedElements_{};
    const vector<double> clusterCutoffs_{};
    SymmetricCE symCE_;
    const vector<vector<vector<int>>> localOrbitsEncoding_;
    const unordered_map<string, double> chemicalPotentialsMap_;
    const unordered_map<string, int> elementCountMap_{};
    const vector<double> ecis_{};
};

#endif // LMC_PRED_INCLUDE_ENERGYPREDICTOR_H_

