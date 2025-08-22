#include "EnergyPredictor.h"

#include "PrintUtility.h"

// This should just take the config and predictorFilename and
// Rest should be read from the param file

EnergyPredictor::EnergyPredictor(
    const ClusterExpansionParameters &ceParams,
    const Config &config) : maxClusterSize_(ceParams.GetMaxClusterSizeCE()),
                            maxBondOrder_(ceParams.GetMaxBondOrderCE()),
                            ecis_(ceParams.GetECIs()),
                            atomicBasis_(
                                ceParams.GetElementSetCE(),
                                ceParams.GetBasisType()),
                            equivalentClustersEncoding_(
                                GetEquivalentClustersEncoding(
                                    config,
                                    maxBondOrder_,
                                    maxClusterSize_))
{
  auto atomVector = config.GetAtomVector();
  set<Element> configElementSet(atomVector.begin(), atomVector.end());

  auto elementSet = ceParams.GetElementSetCE();

  for (auto ele : configElementSet)
  {
    if (elementSet.find(ele) == elementSet.end())
    {
      cerr << "Error in `EnergyPredictor`: element " << ele
           << " is not in the allowed element set!\n";

      cerr << "This cluster expansion is defined only for the following elements: ";
      for (const auto &allowedEle : elementSet)
        cerr << allowedEle << " ";
      cerr << endl;

      exit(1);
    }
  }
}

double EnergyPredictor::ComputeEnergyOfSite(
    const Config &config,
    const size_t &latticeId)
{
  // Can be Stored for each site or may be not
  auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForSite(
      config,
      latticeId,
      maxBondOrder_);

  // Need to include the site as well to be consistent with the encoding
  canonicalSortedLatticeIds.emplace_back(latticeId);

  VectorXd correlationVector = GetCorrelationVector(
      config,
      atomicBasis_,
      canonicalSortedLatticeIds,
      equivalentClustersEncoding_);

  // E = J.Φ_α
  double energyValue = ecis_.dot(correlationVector) / config.GetNumAtoms();

  return energyValue;
}

double EnergyPredictor::ComputeEnergyOfConfig(
    const Config &config)
{
  double totalEnergy = 0;

  for (size_t latticeId = 0; latticeId < config.GetNumLattices(); latticeId++)
  {
    totalEnergy += ComputeEnergyOfSite(config, latticeId);
  }

  return totalEnergy;
}

double EnergyPredictor::GetDeMigration(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  auto canonicalSortedLatticeIdsSite1 = GetCanonicalSortedSitesForSite(
      config,
      latticeIdJumpPair.first,
      maxBondOrder_);

  auto canonicalSortedLatticeIdsSite2 = GetCanonicalSortedSitesForSite(
      config,
      latticeIdJumpPair.second,
      maxBondOrder_);

  cout << "Site1 : ";
  cout << latticeIdJumpPair.first << endl;
  print1DVector(canonicalSortedLatticeIdsSite1);
  cout << "Site2 : ";
  cout << latticeIdJumpPair.second << endl;

  print1DVector(canonicalSortedLatticeIdsSite2);

  return 0;
}

double EnergyPredictor::GetDeSwap(
    Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  // Before Swap

  auto nnLatticeIds = config.GetNeighboringLatticeIdSetOfPair(latticeIdJumpPair, maxBondOrder_);
  nnLatticeIds.insert(latticeIdJumpPair.first);
  nnLatticeIds.insert(latticeIdJumpPair.second);

  double energyBeforeSwap = 0;

  for (const auto latticeId : nnLatticeIds)
  {
    energyBeforeSwap += ComputeEnergyOfSite(config, latticeId);
  }

  config.LatticeJump(latticeIdJumpPair);

  // After Swap
  double energyAfterSwap = 0;

  for (const auto latticeId : nnLatticeIds)
  {
    energyAfterSwap += ComputeEnergyOfSite(config, latticeId);
  }
  
  auto dE = energyAfterSwap - energyBeforeSwap;

  config.LatticeJump(latticeIdJumpPair);

  return dE;
}