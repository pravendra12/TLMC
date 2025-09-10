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
  cout << left << setw(30) << "LVFE ECIs size:" << lvfeECIs_.size() << "\n";

  cout << string(width, '-') << "\n\n";
}

VectorXd LVFEPredictor::GetLocalSiteClusterVector(
    const Config &config,
    const size_t &vacancyLatticeId)
{
  // Will be cached
  auto sortedIdsAroundSite = GetCanonicalSortedSitesForSite(
      config,
      vacancyLatticeId,
      maxBondOrder_);

  VectorXd correlationVector = GetCorrelationVector(
      config,
      atomicBasis_,
      sortedIdsAroundSite,
      encodedOrbitsForSite_);

  return correlationVector;
}

double LVFEPredictor::GetEffectiveVacancyFormationEnergy(
    const Config &config,
    const size_t &vacancyLatticeId)
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
)
{
  // Will be cached
  auto sortedIdsAroundSite = GetCanonicalSortedSitesForSite(
      config,
      vacancyLatticeId,
      maxBondOrder_);

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
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  size_t vacancyLatticeId;
  size_t migratingAtomLatticeId;

  if (config.GetElementOfLattice(latticeIdJumpPair.first) == Element("X"))
  {
    vacancyLatticeId = latticeIdJumpPair.first;
    migratingAtomLatticeId = latticeIdJumpPair.second;
  }
  else
  {
    vacancyLatticeId = latticeIdJumpPair.second;
    migratingAtomLatticeId = latticeIdJumpPair.first;
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

/*
double EnergyPredictor::GetEnergyOfConfigWithVacancy(
    Config &config,
    const size_t &vacancyLatticeId)
{
  // ELVFE : From the model
  double efveValue = GetEffectiveVacancyFormationEnergy(
      config,
      vacancyLatticeId);

  double averageEnergy = 0;
  vector<Element> elementVector = {Element("Mo"), Element("Ta")};

  auto originalElement = config.GetElementOfLattice(vacancyLatticeId);

  for (const auto &element : elementVector)
  {
    config.SetElementOfLattice(vacancyLatticeId, element);
    averageEnergy += ComputeEnergyOfConfig(config);
  }

  config.SetElementOfLattice(vacancyLatticeId, originalElement);

  averageEnergy /= static_cast<double>(elementVector.size());

  // ELVFE = Ev - 1/2 * (EA + EB)

  auto energyWithVacancy = efveValue + averageEnergy;

  return energyWithVacancy;
}

*/

/*
// This function returns the local energy of the site with vacancy
double EnergyPredictor::GetEnergyOfSiteWithVacancy(
    const Config &config,
    const size_t &latticeId)
{
  vector<Element> elementVector = {Element("Mo"), Element("Ta")};

  // This will be upto max cluster cutoff
  auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForSite(
      config,
      latticeId,
      1);
  canonicalSortedLatticeIds.emplace_back(latticeId);

  auto multiElementClusterVector = symCE_.GetMultiElementLocalClusterVector(
      config,
      latticeId,
      elementVector,
      canonicalSortedLatticeIds,
      localOrbitsEncoding_);

  cout << "Mo : ";
  print1DVector(multiElementClusterVector[0]);
  cout << endl;

  cout << "Ta : ";
  print1DVector(multiElementClusterVector[1]);
  cout << endl;

  // In principle the cluster vector will be just opposite in the sign except the
  // first term, will check if it works later

  // ELVFE : From the model
  double efveValue = GetEffectiveVacancyFormationEnergy(
      config,
      latticeId);

  double averageEnergy = 0;

  for (int i = 0; i < multiElementClusterVector.size(); i++)
  {
    averageEnergy += GetTotalEnergy(
        multiElementClusterVector[i]);
  }

  // 1/2 * (EA + EB)
  averageEnergy /= static_cast<double>(multiElementClusterVector.size());

  // Compute the effective local vacancy formation energy (ELVFE):
  // ELVFE = Ev - 1/2 * (EA + EB)
  // where EA and EB are the energies of the same site with elements A and B.

  // Ev
  double energyOfSiteVacancy = efveValue + averageEnergy;

  return energyOfSiteVacancy;
}


*/