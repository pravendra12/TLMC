#include "EnergyPredictor.h"

EnergyPredictor::EnergyPredictor(
    const ClusterExpansionParameters &ceParams,
    const Config &supercellConfig,
    const Config &primitiveConfig) : nAtoms_(supercellConfig.GetNumAtoms()),
                                     allowedElements_(
                                         ceParams.GetAllowedElements()),
                                     clusterCutoffs_(
                                         ceParams.GetClusterCutoffs()),
                                     symCE_(supercellConfig,
                                            primitiveConfig,
                                            allowedElements_,
                                            clusterCutoffs_),
                                     localOrbitsEncoding_(
                                         symCE_.GetLocalOrbitsEncoding()),
                                     chemicalPotentialsMap_(
                                         ceParams.GetChemicalPotentialsMap()),
                                     elementCountMap_(
                                         GetElementCountMap(
                                             supercellConfig)),
                                     ecis_(
                                         ceParams.GetECIs("symCE"))

{
  const int width = 80;
  const int labelWidth = 40;
  const int valueWidth = width - labelWidth - 2;

  // Header
  cout << string(width, '-') << "\n";
  cout << setw((width + 22) / 2) << right << "Energy Predictor Info" << "\n";
  cout << string(width, '-') << "\n";

  // Table
  cout << left << setw(labelWidth) << "Number of atoms:"
            << right << setw(valueWidth) << nAtoms_ << "\n";

  cout << left << setw(labelWidth) << "Allowed elements:"
            << right << setw(valueWidth);
  for (const auto &ele : allowedElements_)
    cout << ele << " ";
  cout << "\n";

  cout << left << setw(labelWidth) << "Cluster cutoffs:"
            << right << setw(valueWidth);
  for (const auto &cut : clusterCutoffs_)
    cout << cut << " ";
  cout << "\n";

  cout << left << setw(labelWidth) << "Local orbits encoding size:"
            << right << setw(valueWidth) << localOrbitsEncoding_.size() << "\n";

  cout << left << setw(labelWidth) << "Chemical potentials:"
            << "\n";
  for (const auto &[ele, mu] : chemicalPotentialsMap_)
  {
    cout << "  " << left << setw(labelWidth - 2) << ele
              << right << setw(valueWidth) << mu << "\n";
  }

  cout << left << setw(labelWidth) << "Sym CE ECIs size:"
            << right << setw(valueWidth) << ecis_.size() << "\n";

  cout << string(width, '-') << "\n\n";

}

// Returns total formation energy
/*
  Formation energy per atom:

      Ef_per_atom = (E_total - (N_A * μ_A + N_B * μ_B)) / (N_A + N_B)

  - E_total: total energy of the supercell
  - N_A, N_B: number of atoms of type A and B
  - μ_A, μ_B: reference chemical potentials of elements A and B

  The cluster expansion model predicts Ef_per_atom,
  but for most calculations (e.g., energy differences or KMC jumps),
  we are interested in the total energy E_total.
*/

double EnergyPredictor::GetTotalFormationEnergy(
    const vector<double> &clusterVector)
{

  double formationEnergy = 0; // per atom

  // print1DVector(ecis_);

  // cout << ecis_.size() << endl;
  // cout << clusterVector.size() << endl;


  for (int i = 0; i < clusterVector.size(); i++)
  {
    formationEnergy += clusterVector[i] * ecis_[i];
  }

  double totalFormationEnergy = nAtoms_ * formationEnergy;

  return totalFormationEnergy;
}

double EnergyPredictor::GetTotalEnergy(
    const vector<double> &clusterVector)
{
  auto totalFormationEnergy = GetTotalFormationEnergy(
      clusterVector);

  double muContribution = 0;
  for (const auto &[element, count] : elementCountMap_)
  {
    // Chemical Potential Value
    auto muVal = chemicalPotentialsMap_.at(element);

    muContribution += count * muVal;
  }

  // Etotal = (N_A + N_B)*Eformation + (N_A * μ_A + N_B * μ_B))
  double totalEnergy = totalFormationEnergy + muContribution;

  return totalEnergy;
}

// This will value will be higher than the formation energy for whole
// supercell due to the constant term contribution to be 1
// But that does not matter as we only dE which can be calculated using
// total formation energy change which is same as dE

double EnergyPredictor::ComputeLocalFormationEnergyOfSite(
    const Config &config,
    const size_t &latticeId)
{
  // Can be Stored for each site or may be not
  auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForSite(
      config,
      latticeId,
      1);

  // Need to include the site as well to be consistent with the encoding
  canonicalSortedLatticeIds.emplace_back(latticeId);

  auto clusterVector = symCE_.GetLocalClusterVector(
      config,
      canonicalSortedLatticeIds,
      localOrbitsEncoding_);

  // E = J.Φ_α
  double energyValue = GetTotalFormationEnergy(clusterVector);

  return energyValue;
}

double EnergyPredictor::ComputeEnergyOfConfig(
    const Config &config)
{
  auto clusterVector = symCE_.GetClusterVectorForConfig(config);

  double totalEnergy = GetTotalEnergy(clusterVector);

  return totalEnergy;
}


double EnergyPredictor::GetDeSwap(
    Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  // Before Swap

  auto energyBeforeSwap = ComputeLocalFormationEnergyOfSite(
                              config,
                              latticeIdJumpPair.first) +
                          ComputeLocalFormationEnergyOfSite(
                              config,
                              latticeIdJumpPair.second);

  config.LatticeJump(latticeIdJumpPair);

  // After Swap
  double energyAfterSwap = ComputeLocalFormationEnergyOfSite(
                               config,
                               latticeIdJumpPair.first) +
                           ComputeLocalFormationEnergyOfSite(
                               config,
                               latticeIdJumpPair.second);

  auto dE = energyAfterSwap - energyBeforeSwap;

  config.LatticeJump(latticeIdJumpPair);

  return dE;
}

unordered_map<string, int> EnergyPredictor::GetElementCountMap(
    const Config &supercellConfig)
{
  unordered_map<string, int> elementCountMap;
  for (auto ele : supercellConfig.GetAtomVector())
  {
    elementCountMap[ele.GetElementString()]++;
  }

  return elementCountMap;
}
