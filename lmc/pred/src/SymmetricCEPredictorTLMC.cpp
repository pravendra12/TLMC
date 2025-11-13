#include "SymmetricCEPredictorTLMC.h"

SymmetricCEPredictorTLMC::SymmetricCEPredictorTLMC(
    const ClusterExpansionParameters &ceParams,
    const TiledSupercell &tiledSupercell,
    const Config &primitiveConfig) : nAtoms_(tiledSupercell.GetTotalNumOfSites()),
                                     allowedElements_(
                                         ceParams.GetAllowedElements()),
                                     clusterCutoffs_(
                                         ceParams.GetClusterCutoffs()),
                                     symCE_(tiledSupercell.GetSmallConfig(),
                                            primitiveConfig,
                                            allowedElements_,
                                            clusterCutoffs_),
                                     localOrbitsEncoding_(
                                         symCE_.GetLocalOrbitsEncoding()),
                                     chemicalPotentialsMap_(
                                         ceParams.GetChemicalPotentialsMap()),
                                     elementCountMap_(
                                         GetElementCountMap(
                                             tiledSupercell.GetSmallConfig())),
                                     ecis_(
                                         ceParams.GetECIs("symCE"))
{
  GetSymmetricallySortedLatticeEncodedMappings(tiledSupercell);
  PrintSymmetricCEPredictorInfo();
}

const vector<string> &SymmetricCEPredictorTLMC::GetAllowedElements() const
{
  return allowedElements_;
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

double SymmetricCEPredictorTLMC::GetTotalFormationEnergy(
    const vector<double> &clusterVector) const
{

  double formationEnergy = 0; // per atom

  for (int i = 0; i < clusterVector.size(); i++)
  {
    formationEnergy += clusterVector[i] * ecis_[i];
  }

  double totalFormationEnergy = nAtoms_ * formationEnergy;

  return totalFormationEnergy;
}

double SymmetricCEPredictorTLMC::GetTotalEnergy(
    const vector<double> &clusterVector) const
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

double SymmetricCEPredictorTLMC::GetTotalEnergy(
    const vector<double> &clusterVector,
    const unordered_map<string, int> &elementCountMap) const
{
  auto totalFormationEnergy = GetTotalFormationEnergy(
      clusterVector);

  double muContribution = 0;
  for (const auto &[element, count] : elementCountMap)
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
// will take <LATTICE_ID, SMALL_CONFIG_IDX>
double SymmetricCEPredictorTLMC::ComputeLocalFormationEnergyOfSite(
    const TiledSupercell &tiledSupercell,
    const LatticeSiteMapping &latticeSiteMapping)
{

  // symmetricLatticeSiteEncodedMappings_ is indexed using latticeId
  // this does not contain latticeId
  auto symmetricSortedLatticeSitesEncoded = symmetricLatticeSiteEncodedMappings_.at(
      latticeSiteMapping.latticeId);

  // This LatticeSiteEncodedMapping for latticeId is not present in
  // symmetricSortedLatticeSitesEncoded hence need to be added and
  //  will lie in the same smallConfig hence assign -1
  symmetricSortedLatticeSitesEncoded.emplace_back(
      LatticeSiteEncodedMapping(
          latticeSiteMapping.latticeId,
          -1));

  auto clusterVector = symCE_.GetLocalClusterVector(
      tiledSupercell,
      latticeSiteMapping.smallConfigId,
      symmetricSortedLatticeSitesEncoded,
      localOrbitsEncoding_);

  cout << clusterVector.size() << endl;

  // E = J.Φ_α
  double energyValue = GetTotalFormationEnergy(clusterVector);

  return energyValue;
}

double SymmetricCEPredictorTLMC::ComputeLocalFormationEnergyForPair(
    const TiledSupercell &tiledSupercell,
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSitePair, // here the elements will be assigned in the same order
    // This is expected to get desired behaviour when
    // latticeSitePairElements.first <-> latticeSitePair.first
    // latticeSitePairElements.second <-> latticeSitePair.second
    const pair<Element, Element> &latticeSitePairElements)
{

  // 1 AT THE END OF VARIABLE CORRESPONDS TO latticeSitePair.first
  // 2 AT THE END OF VARIABLE CORRESPONDS TO latticeSitePair.second

  // symmetricLatticeSiteEncodedMappings_ is indexed using latticeId
  // SITE 1
  size_t smallConfigIdx1 = latticeSitePair.first.smallConfigId;
  auto symmetricSortedLatticeSitesEncoded1 = symmetricLatticeSiteEncodedMappings_.at(
      latticeSitePair.first.latticeId);
  symmetricSortedLatticeSitesEncoded1.emplace_back(
      LatticeSiteEncodedMapping(
          latticeSitePair.first.latticeId,
          -1));

  // SITE 2
  size_t smallConfigIdx2 = latticeSitePair.second.smallConfigId;
  auto symmetricSortedLatticeSitesEncoded2 = symmetricLatticeSiteEncodedMappings_.at(
      latticeSitePair.second.latticeId);
  symmetricSortedLatticeSitesEncoded2.emplace_back(
      LatticeSiteEncodedMapping(
          latticeSitePair.second.latticeId,
          -1));

  // CLUSTER VECTOR
  auto clusterVectorPair = symCE_.GetLocalClusterVectorForPair(
      tiledSupercell,
      latticeSitePair,
      latticeSitePairElements,
      symmetricSortedLatticeSitesEncoded1,
      symmetricSortedLatticeSitesEncoded2,
      localOrbitsEncoding_);

  if (ecis_.size() != clusterVectorPair.first.size() ||
      ecis_.size() != clusterVectorPair.second.size())
  {
    ostringstream msg;
    msg << "Error in `SymmetricCEPredictorTLMC::ComputeLocalFormationEnergyForPair`: "
        << "ECI size (" << ecis_.size() << ") does not match "
        << "clusterVectorPair.first size (" << clusterVectorPair.first.size() << ") or "
        << "clusterVectorPair.second size (" << clusterVectorPair.second.size() << ").";
    throw runtime_error(msg.str());
  }

  // E = J.Φ_α
  double formationEnergy = 0;

  for (int i = 0; i < ecis_.size(); i++)
  {
    // Contribution of first site
    formationEnergy += clusterVectorPair.first[i] * ecis_[i];
    // Contribution of second site
    formationEnergy += clusterVectorPair.second[i] * ecis_[i];
  }

  double totalFormationEnergy = nAtoms_ * formationEnergy;

  return totalFormationEnergy;
}

// Const Object
double SymmetricCEPredictorTLMC::GetDeSwap(
    const TiledSupercell &tiledSupercell,
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeIdJumpPair)
{
  auto element1 = tiledSupercell.GetElementAtSite(latticeIdJumpPair.first);
  auto element2 = tiledSupercell.GetElementAtSite(latticeIdJumpPair.second);

  // Before Swap
  // element1 <-> latticeIdJumpPair.first
  // element2 <-> latticeIdJumpPair.second

  pair<Element, Element> latticeSitePairElements = {element1, element2};

  // Function Swaps it as
  // element1 <-> latticeIdJumpPair.second
  // element2 <-> latticeIdJumpPair.first

  double dEValue = GetDeSwap(
      tiledSupercell,
      latticeIdJumpPair,
      latticeSitePairElements);

  return dEValue;
}

double SymmetricCEPredictorTLMC::GetDeSwap(
    const TiledSupercell &tiledSupercell,
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSitePair,
    const pair<Element, Element> &latticeSitePairElements)
{
  // Before Swap
  // element1 <-> latticeSitePair.first
  // element2 <-> latticeSitePair.second

  auto energyBeforeSwap = ComputeLocalFormationEnergyForPair(
      tiledSupercell,
      latticeSitePair,
      latticeSitePairElements);

  // After Swap
  // element1 will move to latticeSitePair.second
  // element2 will move to latticeSitePair.first

  pair<Element, Element> afterSwapPairElements = {
      latticeSitePairElements.second,
      latticeSitePairElements.first};

  double energyAfterSwap = ComputeLocalFormationEnergyForPair(
      tiledSupercell,
      latticeSitePair,
      afterSwapPairElements);

  auto dE = energyAfterSwap - energyBeforeSwap;

  return dE;
}

unordered_map<string, int> SymmetricCEPredictorTLMC::GetElementCountMap(
    const Config &supercellConfig)
{
  unordered_map<string, int> elementCountMap;
  for (auto ele : supercellConfig.GetAtomVector())
  {
    elementCountMap[ele.GetElementString()]++;
  }

  return elementCountMap;
}

void SymmetricCEPredictorTLMC::GetSymmetricallySortedLatticeEncodedMappings(
    const TiledSupercell &tiledSupercell)
{
  // Tiled suprecell neighbour list must also be updated as well

  // Expecting supercellConfig to be updated in the following way
  // supercellConfig.UpdateNeighbours({maxClusterCutoff})
  const size_t maxBondOrder = 1;
  const Config &smallConfig = tiledSupercell.GetSmallConfig();

  // this is reference to avoid copy
  // double maxClusterCutoff = *std::max_element(clusterCutoffs_.begin(), clusterCutoffs_.end());
  // smallConfig.UpdateNeighborList({maxClusterCutoff});

  // Number of lattice in the small config
  const size_t numLattices = smallConfig.GetNumLattices();

  symmetricLatticeSiteEncodedMappings_.reserve(numLattices);

  size_t numSortedIds = smallConfig.GetNeighborLatticeIdVectorOfLattice(0, 1).size();

  vector<LatticeSiteEncodedMapping> canonicalSortedLatticeAndConfigIdVector;
  canonicalSortedLatticeAndConfigIdVector.reserve(numSortedIds);

  unordered_map<size_t, size_t> latticeIdToSortedIndexMap;
  latticeIdToSortedIndexMap.reserve(numSortedIds);

  // Iterate over all the sites and get store the sorted lattice Ids
  for (size_t latticeId = 0; latticeId < numLattices; latticeId++)
  {
    auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForSite(
        smallConfig,
        latticeId,
        maxBondOrder);

    // this latticeId will not be there in the neighbourlist as one need to merge it with a
    // smallConfigIdx
    // canonicalSortedLatticeIds.emplace_back(latticeId);

    // Clear to avoid any unwanted issues
    latticeIdToSortedIndexMap.clear();
    for (size_t i = 0; i < numSortedIds; i++)
    {
      latticeIdToSortedIndexMap[canonicalSortedLatticeIds[i]] = i;
    }

    // now get the neighbourlist for the same lattice Id from the tiledSupercell
    // and build the canonicalSortedLatticeAndConfigIdVector
    // <nnLatticeId, cubeIdx>
    canonicalSortedLatticeAndConfigIdVector.clear();
    canonicalSortedLatticeAndConfigIdVector = tiledSupercell.GetNeighborLatticeIdVectorOfLattice(
        latticeId,
        maxBondOrder);

    // Sort neighboursInTiledSupercell based on the order in `order`
    std::sort(canonicalSortedLatticeAndConfigIdVector.begin(), canonicalSortedLatticeAndConfigIdVector.end(),
              [&](const auto &a, const auto &b)
              {
                return latticeIdToSortedIndexMap.at(a.latticeId) < latticeIdToSortedIndexMap.at(b.latticeId);
              });

    symmetricLatticeSiteEncodedMappings_.emplace_back(move(canonicalSortedLatticeAndConfigIdVector));
  }
}

void SymmetricCEPredictorTLMC::PrintSymmetricCEPredictorInfo() const
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

  cout << left << setw(labelWidth) << "Cluster vector size (with constant term):"
       << right << setw(valueWidth) << (localOrbitsEncoding_.size() + 1) << "\n";

  cout << left << setw(labelWidth) << "Chemical potentials:"
       << "\n";
  for (const auto &[ele, mu] : chemicalPotentialsMap_)
  {
    cout << "  " << left << setw(labelWidth - 2) << ele
         << right << setw(valueWidth) << mu << "\n";
  }
  cout << left << setw(labelWidth) << "Size of precomputed sortedIds vector:"
       << symmetricLatticeSiteEncodedMappings_.size() << "\n";

  cout << left << setw(labelWidth) << "Sym CE ECIs size:"
       << right << setw(valueWidth) << ecis_.size() << "\n";

  cout << string(width, '-') << "\n\n";
}