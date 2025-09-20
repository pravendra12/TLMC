#include "LVFEPredictorTLMC.h"

LVFEPredictorTLMC::LVFEPredictorTLMC(
    const ClusterExpansionParameters &ceParams,
    const TiledSupercell &tiledSupercell) : maxBondOrder_(ceParams.GetMaxBondOrder("lvfe")),
                                            maxClusterSize_(
                                                ceParams.GetMaxClusterSize("lvfe")),
                                            atomicBasis_(
                                                ceParams.GetElementSet("lvfe"),
                                                ceParams.GetBasisType()),
                                            encodedOrbitsForSite_(
                                                GetLocalEncodedOrbitsForSite(
                                                    tiledSupercell.GetSmallConfig(),
                                                    maxBondOrder_,
                                                    maxClusterSize_,
                                                    true)),
                                            lvfeECIs_(ceParams.GetECIs("lvfe"))

{
  // Precompute the sorted lattice Ids
  GetSymmetricallySortedLatticeIdsVectorMap(tiledSupercell);
  // Print the info
  PrintLVFEPredictor();
}

VectorXd LVFEPredictorTLMC::GetLocalSiteClusterVector(
    const TiledSupercell &tiledSupercell,
    const LatticeSiteMapping &vacancyLatticeSite) const
{

  // No need to include the vacancyLatticeSite
  auto sortedIdsAroundSite = symmetricLatticeSiteEncodedMappings_[vacancyLatticeSite.latticeId];

  VectorXd correlationVector = GetCorrelationVector(
      tiledSupercell,
      vacancyLatticeSite,
      atomicBasis_,
      sortedIdsAroundSite,
      encodedOrbitsForSite_);

  return correlationVector;
}

double LVFEPredictorTLMC::GetEffectiveVacancyFormationEnergy(
    const TiledSupercell &tiledSupercell,
    const LatticeSiteMapping &vacancyLatticeSite) const
{
  VectorXd correlationVector = GetLocalSiteClusterVector(
      tiledSupercell,
      vacancyLatticeSite);

  // local vacancy formation energy
  double lvfeValue = 0;

  for (int i = 0; i < lvfeECIs_.size(); i++)
  {
    lvfeValue += lvfeECIs_[i] * correlationVector(i);
  }

  return lvfeValue;
}

double LVFEPredictorTLMC::GetEffectiveVacancyFormationEnergy(
    const TiledSupercell &tiledSupercell,
    const LatticeSiteMapping &vacancyLatticeSite,
    const LatticeSiteMapping &targetLatticeSite, // LatticeId to which the element will be assigned
    const Element &elementToAssign               // Element which need to be assigned to the lattice Id
) const
{
  auto sortedIdsAroundSite = symmetricLatticeSiteEncodedMappings_[vacancyLatticeSite.latticeId];

  // Will return the correlation vector with element assigned to a targetLatticeSite
  VectorXd correlationVector = GetCorrelationVector(
      tiledSupercell,
      vacancyLatticeSite,
      atomicBasis_,
      targetLatticeSite,
      elementToAssign,
      sortedIdsAroundSite,
      encodedOrbitsForSite_);

  // local vacancy formation energy
  double lvfeValue = 0;

  for (int i = 0; i < lvfeECIs_.size(); i++)
  {
    lvfeValue += lvfeECIs_[i] * correlationVector(i);
  }

  return lvfeValue;
}

double LVFEPredictorTLMC::GetDeForVacancyMigration(
    const TiledSupercell &tiledSupercell,
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair) const
{
  LatticeSiteMapping vacancyLatticeSite;
  LatticeSiteMapping migratingAtomLatticeSite;

  bool firstIsVacancy = (tiledSupercell.GetElementAtSite(latticeSiteJumpPair.first) == Element("X"));
  bool secondIsVacancy = (tiledSupercell.GetElementAtSite(latticeSiteJumpPair.second) == Element("X"));

  if (firstIsVacancy && !secondIsVacancy)
  {
    vacancyLatticeSite = latticeSiteJumpPair.first;
    migratingAtomLatticeSite = latticeSiteJumpPair.second;
  }
  else if (!firstIsVacancy && secondIsVacancy)
  {
    vacancyLatticeSite = latticeSiteJumpPair.second;
    migratingAtomLatticeSite = latticeSiteJumpPair.first;
  }
  else
  {
    throw std::runtime_error("Error in `LVFEPredictorTLMC::GetDeForVacancyMigration`: One of the lattice sites must be a vacancy.");
  }

  auto migratingElement = tiledSupercell.GetElementAtSite(migratingAtomLatticeSite);

  // Energy Before Vacancy Jump
  // Ev before the jump

  // Effective vacancy formation energy
  auto lvfeValueBeforeJump = GetEffectiveVacancyFormationEnergy(
      tiledSupercell,
      vacancyLatticeSite);

  // Energy After Vacancy Jump
  // migratingElement will migrate to vacancyLatticeSite
  // migratingAtomLatticeId will be the new vacancyLatticeSite

  auto lvfeValueAfterJump = GetEffectiveVacancyFormationEnergy(
      tiledSupercell,
      migratingAtomLatticeSite, // newvacancyLatticeSite
      vacancyLatticeSite,       // migrating element will be assigned to vacancyLatticeSite
      migratingElement);

  // ELVFE = Ev - 1/2 * (EA + EB)
  // For a site EA + EB with 1 -1 will lead to only 2*Jo
  // Same for the other site EA + EB with 1 -1 will lead to only 2*Jo
  // Leading to cancelation of Jo terms

  double dEValue = lvfeValueAfterJump - lvfeValueBeforeJump;

  return dEValue;
}

void LVFEPredictorTLMC::GetSymmetricallySortedLatticeIdsVectorMap(
    const TiledSupercell &tiledSupercell)
{
  const Config &smallConfig = tiledSupercell.GetSmallConfig();

  // Number of lattice in the small config
  const size_t numLattices = smallConfig.GetNumLattices();

  symmetricLatticeSiteEncodedMappings_.reserve(numLattices);

  auto numSortedSites = smallConfig.GetNeighborLatticeIdsUpToOrder(0, maxBondOrder_).size();
  vector<LatticeSiteEncodedMapping> sortedLatticeSitesEncoded;
  sortedLatticeSitesEncoded.reserve(numSortedSites);

  unordered_map<size_t, size_t> latticeIdToSortedIndexMap;
  latticeIdToSortedIndexMap.reserve(numSortedSites);

  // Iterate over all the sites and get store the sorted lattice Ids
  for (size_t latticeId = 0; latticeId < numLattices; latticeId++)
  {
    auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForSite(
        smallConfig,
        latticeId,
        maxBondOrder_);

    // Need to be cleared as the sorted lattice sites will be sorted in a unique order
    // for a give site
    latticeIdToSortedIndexMap.clear();
    for (size_t i = 0; i < numSortedSites; i++)
    {
      latticeIdToSortedIndexMap[canonicalSortedLatticeIds[i]] = i;
    }

    // now get the neighbourlist for the same lattice Id from the tiledSupercell
    // and build the canonicalSortedLatticeAndConfigIdVector
    // <nnLatticeId, cubeIdx>
    sortedLatticeSitesEncoded.clear();
    sortedLatticeSitesEncoded = tiledSupercell.GetNeighborLatticeIdsUpToOrder(
        latticeId,
        maxBondOrder_);

    std::sort(sortedLatticeSitesEncoded.begin(), sortedLatticeSitesEncoded.end(),
              [&](const auto &a, const auto &b)
              {
                return latticeIdToSortedIndexMap.at(a.latticeId) < latticeIdToSortedIndexMap.at(b.latticeId);
              });

    symmetricLatticeSiteEncodedMappings_.emplace_back(move(sortedLatticeSitesEncoded));
  }
}

void LVFEPredictorTLMC::PrintLVFEPredictor() const
{
  const int width = 80;
  const int labelWidth = 30;
  const int valueWidth = width - labelWidth - 2; // 2 for spacing

  // Header
  cout << string(width, '-') << "\n";
  cout << setw((width + 22) / 2) << right << "LVFE Predictor Info" << "\n";
  cout << string(width, '-') << "\n";

  cout << left << setw(30) << "Max bond order:" << maxBondOrder_ << "\n";
  cout << left << setw(30) << "Max cluster size:" << maxClusterSize_ << "\n";
  cout << left << setw(30) << "Encoded orbits size:" << encodedOrbitsForSite_.size() << "\n";
  cout << left << setw(30) << "Size of precomputed sortedIds vector:"
       << symmetricLatticeSiteEncodedMappings_.size() << "\n";
  cout << left << setw(30) << "LVFE ECIs size:" << lvfeECIs_.size() << "\n";

  cout << string(width, '-') << "\n\n";
}