// main for Tiled KMC

#include "Config.h"
#include "Cube.h"
#include "TiledSupercell.h"
#include <filesystem>
#include "TiledSupercell.h"
#include "Config.h"
namespace fs = filesystem;

int main()
{
  size_t smallConfigSize = 10;
  double latticeParam = 3.2;
  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");

  Cube cubeObj(3);

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

  auto largeConfig = Config::GenerateSupercell(30, 3.2, "Mo", "BCC");
  
  for (int i = 0; i < 10000; i++)
  {
    largeConfig.SetElementOfLattice(i, Element("Ta"));
  }
  largeConfig.SetElementOfLattice(0, Element("X"));

  cout << "Done ReadCfg" << endl;
  
  tiledSupercell.UpdateAtomVector(
      largeConfig);

  cout << "Done UpdateAtomVector" << endl;

  Config::WriteConfig(
      "testingAtomIndicesReadWrite/originalConfig.cfg.gz", largeConfig);

  cout << "Done WriteConfig" << endl;

  auto largeCfgFromTS = tiledSupercell.MakeSupercell();

  Config::WriteConfig(
      "testingAtomIndicesReadWrite/configFromMakeSupercell.cfg.gz", largeConfig);

  cout << "Done WriteConfig" << endl;

  tiledSupercell.WriteAtomicIndicesToFile("testingAtomIndicesReadWrite/atomIndices.bin.gz");

  cout << "Done WriteAtomicIndicesToFile" << endl;

  vector<uint64_t> atomicIndicesVector = TiledSupercell::ReadAtomicIndicesFromFile("testingAtomIndicesReadWrite/atomIndices.bin.gz");

  tiledSupercell.UpdateAtomVector(atomicIndicesVector);

  auto largeCfgFromAtomicIndices = tiledSupercell.MakeSupercell();

  Config::WriteConfig(
      "testingAtomIndicesReadWrite/configFromAtomicIndices.cfg.gz", largeConfig);

  const string pathToAtomVector = "//media/sf_Phd/MoTa/TLMC/kmc_1000K/initalConfigFromCMC";
  const string pathToOutput = "//media/sf_Phd/MoTa/TLMC/kmc_1000K/initalConfigFromCMC/configs";

  exit(1);

  for (const auto &entry : fs::directory_iterator(pathToAtomVector))
  {
    if (!entry.is_regular_file())
      continue;

    auto filePath = entry.path();
    if (filePath.extension() != ".txt")
      continue; // only .txt files

    string baseName = filePath.stem().string(); // e.g. "0" from "0.txt"

    // check if basename is an integer
    bool isInt = true;
    try
    {
      size_t pos;
      stoi(baseName, &pos);
      if (pos != baseName.size())
        isInt = false; // extra chars after number
    }
    catch (...)
    {
      isInt = false;
    }
    if (!isInt)
      continue;

    // read atomVector from file
    auto atomVector = TiledSupercell::ReadAtomVectorInfoFromFile(filePath.string());

    // build large config from atomVector
    auto largeCfgFromAtomVector = TiledSupercell::MakeSupercellFromAtomInfo(
        smallCfg, cubeObj, atomVector);

    // output filename
    string outFile = pathToOutput + "/" + baseName + ".cfg.gz";

    Config::WriteConfig(outFile, largeCfgFromAtomVector);
    cout << "Wrote " << outFile << endl;
  }
}

/*
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
*/