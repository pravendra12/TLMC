// main for Tiled KMC

#include "Config.h"
#include "Cube.h"
#include "TiledSupercell.h"
#include "EnergyPredictor.h"
#include "SymmetricCEPredictorTLMC.h"

#include "LVFEPredictor.h"
#include "LVFEPredictorTLMC.h"
#include "KRAPredictor.h"
#include "KRAPredictorTLMC.h"
#include "VacancyMigrationPredictorTLMC.h"
#include "EnergyPredictorTLMC.h"

int main()
{
  size_t smallConfigSize = 10;
  double latticeParam = 3.2;
  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");

  smallCfg.UpdateNeighborList({3});
  smallCfg.SetElementOfLattice(1, Element("X"));
  smallCfg.SetElementOfLattice(
      smallCfg.GetNeighborLatticeIdVectorOfLattice(1, 1)[2], Element("Ta"));
  smallCfg.SetElementOfLattice(
      smallCfg.GetNeighborLatticeIdVectorOfLattice(1, 1)[3], Element("Ta"));

  auto primCfg = Config::GenerateSupercell(
      2,
      3.2,
      "Mo",
      "BCC");

  size_t configIdx = 0;

  size_t vacancyId = smallCfg.GetVacancyLatticeId();
  auto vacancySiteId = LatticeSiteMapping(
      vacancyId,
      configIdx);

  auto nnLatticeId = smallCfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[0];

  auto nnLatticeSiteId = LatticeSiteMapping(
      nnLatticeId,
      configIdx);

  cout << "VacancyId : " << vacancyId << endl;
  cout << "nnLatticeId : " << nnLatticeId << endl;

  // Now to make a 50x50x50 supercell with 10x10x10 small suprecells
  // One would need to stack 5 of these small cfg in a cube
  // total 5x5x5 cube = 125 cubes
  // Total Atoms = 125 * 2000
  Cube cubeObj(2);

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

  ClusterExpansionParameters ceParams(
      "/home/pravendra3/Documents/TLMC/script/coefficientFile_MoTa_V3.1.json");

  cout << tiledSupercell.GetElementAtSite({0, 0}).GetElementString() << endl;
  cout << tiledSupercell.GetElementAtSite({1, 0}).GetElementString() << endl;

  smallCfg.UpdateNeighborList({13}); // update using max cluster cutoff

  tiledSupercell.UpdateNeighbourLists(1);
  SymmetricCEPredictorTLMC symCEPredictorTLMC(
      ceParams,
      tiledSupercell,
      primCfg);

  smallCfg.UpdateNeighborList({3, 4, 5});
  tiledSupercell.UpdateNeighbourLists(3);

  LVFEPredictorTLMC lvfePredictorTLMC(
      ceParams,
      tiledSupercell);

  KRAPredictorTLMC kraPredictorTLMC(
      ceParams,
      tiledSupercell);

  EnergyPredictorTLMC energyPredictor(
      symCEPredictorTLMC,
      lvfePredictorTLMC);

  VacancyMigrationPredictorTLMC vacMigrationPredictor(
      kraPredictorTLMC,
      energyPredictor);

  auto dEValue = energyPredictor.GetEnergyChange(
      tiledSupercell,
      make_pair(
          vacancySiteId,
          nnLatticeSiteId));

  auto barrierAndDe = vacMigrationPredictor.GetBarrierAndDeltaE(
      tiledSupercell,
      make_pair(
          vacancySiteId,
          nnLatticeSiteId));

  cout << "Energy Change: " << dEValue << endl;
  cout << "Barrier : " << barrierAndDe.first << endl;
  cout << "dE: " << barrierAndDe.second << endl;

  return 0;
}