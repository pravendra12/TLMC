#include "ComputeVFE.h"

map<Element, double> ComputeVFE(
    Config &config,
    const unordered_map<string, double> &chemicalPotentialMap,
    SymmetricCEPredictor &symCEPredictor,
    LVFEPredictor &lvfePredictor)
{
  const size_t centralLatticeId = config.GetCentralAtomLatticeId();
  const auto centralAtomElement = config.GetElementOfLattice(centralLatticeId);

  // Compute the LVFE
  double energyLVFE = lvfePredictor.GetEffectiveVacancyFormationEnergy(config, centralLatticeId);

  double avgEnergy = 0.0;
  map<Element, double> energyWithElementMap;

  for (const auto &elementEntry : chemicalPotentialMap)
  {
    const auto element = Element(elementEntry.first);

    config.SetElementOfLattice(centralLatticeId, element);
    double energyOfConfig = symCEPredictor.ComputeEnergyOfConfig(config);

    energyWithElementMap[element] = energyOfConfig;
    avgEnergy += energyOfConfig;

    //cout << "energyOfConfig : " << energyOfConfig << endl;

    // Restore the original element
    config.SetElementOfLattice(centralLatticeId, centralAtomElement);
  }

  avgEnergy /= double(chemicalPotentialMap.size());

  //cout << "ELVFE : " << energyLVFE << endl;

  // Elvfe = Ev - avgEnergy
  double energyWithVacancy = energyLVFE + avgEnergy;

  cout << "energyWithVacancy : " << energyWithVacancy << endl;


  map<Element, double> vfeMap;
  for (const auto &elementEntry : energyWithElementMap)
  {
    vfeMap[elementEntry.first] = energyWithVacancy - elementEntry.second + chemicalPotentialMap.at(elementEntry.first.GetElementString());
  }

  return vfeMap;
}
