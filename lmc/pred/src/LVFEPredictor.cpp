#include "LVFEPredictor.h"

LVFEPredictor::LVFEPredictor(
    const ClusterExpansionParameters &ceParams,
    const Config &config) : maxBondOrder_(ceParams.GetMaxBondOrder("lvfe")),
                            maxClusterSize_(
                                ceParams.GetMaxClusterSize("lvfe")),
                            atomicBasis_(
                                ceParams.GetElementSet("lvfe"),
                                ceParams.GetBasisType()),
                            encodedOrbitsForSite_(
                                GetLocalEncodedOrbitsForSite(
                                    config,
                                    maxBondOrder_,
                                    maxClusterSize_,
                                    true)),
                            lvfeECIs_(ceParams.GetECIs("lvfe"))

{
  // Precompute the sorted lattice Ids
  GetSymmetricallySortedLatticeIdsVectorMap(config);
  // Print the info
  PrintLVFEPredictor();
}

VectorXd LVFEPredictor::GetLocalSiteClusterVector(
    const Config &config,
    const size_t &vacancyLatticeId) const
{
  // Will be cached
  /*
  auto sortedIdsAroundSite = GetCanonicalSortedSitesForSite(
      config,
      vacancyLatticeId,
      maxBondOrder_);
  */

  // No need to include the vacancyLatticeId
  auto sortedIdsAroundSite = symmetricallySortedLatticeIdsVectorMap_[vacancyLatticeId];

  VectorXd correlationVector = GetCorrelationVector(
      config,
      atomicBasis_,
      sortedIdsAroundSite,
      encodedOrbitsForSite_);

  return correlationVector;
}

double LVFEPredictor::GetEffectiveVacancyFormationEnergy(
    const Config &config,
    const size_t &vacancyLatticeId) const
{
  VectorXd correlationVector = GetLocalSiteClusterVector(
      config,
      vacancyLatticeId);

  // local vacancy formation energy
  double lvfeValue = 0;

  for (int i = 0; i < lvfeECIs_.size(); i++)
  {
    lvfeValue += lvfeECIs_[i] * correlationVector(i);
  }

  return lvfeValue;
}

double LVFEPredictor::GetEffectiveVacancyFormationEnergy(
    const Config &config,
    const size_t &vacancyLatticeId,
    const size_t &targetLatticeId, // LatticeId to which the element will be assigned
    const Element &elementToAssign // Element which need to be assigned to the lattice Id
) const
{
  // Will be cached
  /*
  auto sortedIdsAroundSite = GetCanonicalSortedSitesForSite(
      config,
      vacancyLatticeId,
      maxBondOrder_);
  */

  auto sortedIdsAroundSite = symmetricallySortedLatticeIdsVectorMap_[vacancyLatticeId];

  // Will return the correlation vector with element assigned to a targetLatticeId
  VectorXd correlationVector = GetCorrelationVector(
      config,
      atomicBasis_,
      targetLatticeId,
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

double LVFEPredictor::GetDeForVacancyMigration(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair) const
{
  size_t vacancyLatticeId;
  size_t migratingAtomLatticeId;

  bool firstIsVacancy = (config.GetElementOfLattice(latticeIdJumpPair.first) == Element("X"));
  bool secondIsVacancy = (config.GetElementOfLattice(latticeIdJumpPair.second) == Element("X"));

  if (firstIsVacancy && !secondIsVacancy)
  {
    vacancyLatticeId = latticeIdJumpPair.first;
    migratingAtomLatticeId = latticeIdJumpPair.second;
  }
  else if (!firstIsVacancy && secondIsVacancy)
  {
    vacancyLatticeId = latticeIdJumpPair.second;
    migratingAtomLatticeId = latticeIdJumpPair.first;
  }
  else
  {
    throw std::runtime_error("Error in `LVFEPredictor::GetDeForVacancyMigration`: One of the lattice sites must be a vacancy.");
  }

  auto migratingElement = config.GetElementOfLattice(migratingAtomLatticeId);

  // Energy Before Vacancy Jump
  // Ev before the jump

  // Effective vacancy formation energy
  auto lvfeValueBeforeJump = GetEffectiveVacancyFormationEnergy(
      config,
      vacancyLatticeId);

  // Energy After Vacancy Jump
  // migratingElement will migrate to vacancyLatticeId
  // migratingAtomLatticeId will be the new vacancyLatticeId

  auto lvfeValueAfterJump = GetEffectiveVacancyFormationEnergy(
      config,
      migratingAtomLatticeId, // newVacancyLatticeId
      vacancyLatticeId,       // migrating element will be assigned to vacancyLatticeId
      migratingElement);

  // ELVFE = Ev - 1/2 * (EA + EB)
  // For a site EA + EB with 1 -1 will lead to only 2*Jo
  // Same for the other site EA + EB with 1 -1 will lead to only 2*Jo
  // Leading to cancelation of Jo terms

  double dEValue = lvfeValueAfterJump - lvfeValueBeforeJump;

  return dEValue;
}

double LVFEPredictor::GetDeSwap(
    Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair) const
{
  size_t vacancyLatticeId;
  size_t migratingAtomLatticeId;

  bool firstIsVacancy = (config.GetElementOfLattice(latticeIdJumpPair.first) == Element("X"));
  bool secondIsVacancy = (config.GetElementOfLattice(latticeIdJumpPair.second) == Element("X"));

  if (firstIsVacancy && !secondIsVacancy)
  {
    vacancyLatticeId = latticeIdJumpPair.first;
    migratingAtomLatticeId = latticeIdJumpPair.second;
  }
  else if (!firstIsVacancy && secondIsVacancy)
  {
    vacancyLatticeId = latticeIdJumpPair.second;
    migratingAtomLatticeId = latticeIdJumpPair.first;
  }
  else
  {
    throw std::runtime_error("Error in `LVFEPredictor::GetDeSwap`: One of the lattice sites must be a vacancy.");
  }
  auto migratingElement = config.GetElementOfLattice(migratingAtomLatticeId);

  // Energy Before Vacancy Jump
  // Ev before the jump

  // Effective vacancy formation energy
  auto lvfeValueBeforeJump = GetEffectiveVacancyFormationEnergy(
      config,
      vacancyLatticeId);

  config.LatticeJump(latticeIdJumpPair);

  // Energy After Vacancy Jump
  // migratingElement will migrate to vacancyLatticeId
  // migratingAtomLatticeId will be the new vacancyLatticeId

  auto lvfeValueAfterJump = GetEffectiveVacancyFormationEnergy(
      config,
      migratingAtomLatticeId); // vacancy have migrated

  // ELVFE = Ev - 1/2 * (EA + EB)
  // For a site EA + EB with 1 -1 will lead to only 2*Jo
  // Same for the other site EA + EB with 1 -1 will lead to only 2*Jo
  // Leading to cancelation of Jo terms

  double dEValue = lvfeValueAfterJump - lvfeValueBeforeJump;

  config.LatticeJump(latticeIdJumpPair);

  return dEValue;
}

void LVFEPredictor::GetSymmetricallySortedLatticeIdsVectorMap(
    const Config &config) 
{
  const size_t numLattices = config.GetNumLattices();

  symmetricallySortedLatticeIdsVectorMap_.reserve(numLattices);

  auto numSortedSites = config.GetNeighborLatticeIdsUpToOrder(0, maxBondOrder_).size();
  vector<size_t> canonicalSortedLatticeIds;
  canonicalSortedLatticeIds.reserve(numSortedSites);

  // Iterate over all the sites and get store the sorted lattice Ids
  for (size_t latticeId = 0; latticeId < numLattices; latticeId++)
  {
    canonicalSortedLatticeIds = GetCanonicalSortedSitesForSite(
        config,
        latticeId,
        maxBondOrder_);

    symmetricallySortedLatticeIdsVectorMap_.emplace_back(move(canonicalSortedLatticeIds));
  }
}

void LVFEPredictor::PrintLVFEPredictor() const
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
                           << symmetricallySortedLatticeIdsVectorMap_.size() << "\n";
  cout << left << setw(30) << "LVFE ECIs size:" << lvfeECIs_.size() << "\n";

  cout << string(width, '-') << "\n\n";
}