#include "ComputeVFE.h"

map<Element, double> ComputeVFE(
    Config &config,
    const unordered_map<string, double> &chemicalPotentialMap,
    SymmetricCEPredictor &symCEPredictor,
    LVFEPredictor *lvfePredictor)
{
  const size_t centralLatticeId = config.GetCentralAtomLatticeId();
  const auto centralAtomElement = config.GetElementOfLattice(centralLatticeId);

  map<Element, double> energyWithElementMap;
  double energyWithVacancy = 0;
  map<Element, double> vfeMap;

  if (!lvfePredictor)
  {

    // Vacancy formation energy using CE *formation energy per site*.
    //
    // CE was trained on formation energy per site:
    //   e_form(cfg) = (E_abs(cfg) - sum_i N_i * u_i) / N
    // so:
    //   E_abs(cfg) = N * e_form(cfg) + sum_i N_i * u_i
    //
    // If we create a vacancy by replacing one A atom with X (same total sites N):
    //   cfg_A : site = A
    //   cfg_X : site = X   (vacancy species)
    //
    // Then the vacancy formation energy is:
    //   E_vf(A) = E_abs(cfg_X) - E_abs(cfg_A) + mu_A
    //
    // If we take the reservoir chemical potential as the *same pure reference* used in e_form
    // (i.e., mu_A = u_A, typically 0 K pure-element DFT), the u-terms cancel and:
    //
    //   E_vf(A) = N * ( e_form(cfg_X) - e_form(cfg_A) )
    //
    // Note: e_form is per-site, so multiply by N to get total energy (eV).

    for (const auto &elementEntry : chemicalPotentialMap)
    {
      const auto element = Element(elementEntry.first);

      config.SetElementOfLattice(centralLatticeId, element);

      // either compute the energy of entire config or just the local config will give nearly same results
      auto localFormationEnergy = symCEPredictor.ComputeLocalFormationEnergyOfSite(config, centralLatticeId);
      energyWithElementMap[element] = localFormationEnergy;
    }

    // energy with vacancy
    config.SetElementOfLattice(centralLatticeId, Element("X"));
    auto localFormationEnergyWithX = symCEPredictor.ComputeLocalFormationEnergyOfSite(config, centralLatticeId);

    for (const auto &elementEntry : energyWithElementMap)
    {
      // Check the comment above
      vfeMap[elementEntry.first] = localFormationEnergyWithX - elementEntry.second;
    }
  }

  else
  {
    // Compute the LVFE
    double energyLVFE = lvfePredictor->GetEffectiveVacancyFormationEnergy(config, centralLatticeId);

    double avgEnergy = 0.0;

    for (const auto &elementEntry : chemicalPotentialMap)
    {
      const auto element = Element(elementEntry.first);

      config.SetElementOfLattice(centralLatticeId, element);
      double energyOfConfig = symCEPredictor.ComputeEnergyOfConfig(config);

      energyWithElementMap[element] = energyOfConfig;
      avgEnergy += energyOfConfig;

      // Restore the original element
      config.SetElementOfLattice(centralLatticeId, centralAtomElement);
    }

    avgEnergy /= double(chemicalPotentialMap.size());

    // cout << "ELVFE : " << energyLVFE << endl;

    // Elvfe = Ev - avgEnergy
    energyWithVacancy = energyLVFE + avgEnergy;

    // cout << "energyWithVacancy : " << energyWithVacancy << endl;

    for (const auto &elementEntry : energyWithElementMap)
    {
      vfeMap[elementEntry.first] = energyWithVacancy - elementEntry.second + chemicalPotentialMap.at(elementEntry.first.GetElementString());
    }
  }

  config.SetElementOfLattice(centralLatticeId, centralAtomElement);
  return vfeMap;
}
