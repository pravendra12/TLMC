#include "Home.h"
#include "Config.h"
#include "TiledSupercell.h"
#include "Cube.h"
#include <random>

int main()
{
  // CONVERTING TO CONFIG FROM ATOM INDICIES.bin.gz
  size_t smallConfigSize = 10;
  double latticeParam = 3.2;
  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");

  Cube cubeObj(5);

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

  // Converting the atomicVector.txt to configs

  // const string pathToAtomVector = "/media/sf_Phd/codeTLMC/benchmark/kmc1200K/TLMC";
  const string pathToAtomVector = "/media/sf_Phd/codeTLMC/benchmark/cmc1700K/TLMC";

  const string pathToOutput = "/media/sf_Phd/codeTLMC/benchmark/cmc1700K/TLMC/configs";

  for (const auto &entry : fs::directory_iterator(pathToAtomVector))
  {
    // cout << entry.path().string() << endl;

    auto filePath = entry.path();

    // process only .bin.gz files
    if (filePath.extension() != ".gz" || filePath.stem().extension() != ".bin")
      continue;

    string baseName = filePath.stem().stem().string(); // "0" from "0.bin.gz"


    // read atomVector from file
    // auto atomVector = TiledSupercell::ReadAtomVectorInfoFromFile(filePath.string());
    auto atomIndicesVector = TiledSupercell::ReadAtomicIndicesFromFile(filePath.string());

    vector<Element> atomVector;
    atomVector.reserve(atomIndicesVector.size());

    for (const auto atomIdx : atomIndicesVector)
    {
      atomVector.emplace_back(Element(atomIdx));
    }

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
// CONVERSION FROM atomVector.txt file
int main()
{
  size_t smallConfigSize = 10;
  double latticeParam = 3.2;
  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");

  Cube cubeObj(5);

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

  // Converting the atomicVector.txt to configs

  const string pathToAtomVector = "/media/sf_Phd/codeTLMC/benchmark/kmc600K/TLMC";
  const string pathToOutput = "/media/sf_Phd/codeTLMC/benchmark/kmc600K/TLMC/configs";


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
  */

/*
//// MAKING RANDOM CONFIG FOR TILED LATTICE MONTE CARLO SIMULATION
int main()
{

  size_t smallConfigSize = 10;
  double latticeParam = 3.2;
  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");

  Cube cubeObj(30);

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

  // make a atom vector

  // Number of sites
  size_t numSites = static_cast<size_t>(54e6);

  // Prepare vector
  std::vector<uint64_t> atomIndexVector(numSites);

  uint64_t moIndex = Element("Mo").GetAtomicIndex();
  uint64_t taIndex = Element("Ta").GetAtomicIndex();

  // First site is always Mo
  atomIndexVector[0] = uint64_t(0);

  // Fill remaining sites
  size_t half = (numSites - 1) / 2; // minus first site
  for (size_t i = 1; i <= half; i++)
    atomIndexVector[i] = moIndex;
  for (size_t i = half + 1; i < numSites; i++)
    atomIndexVector[i] = taIndex;

  // Shuffle from index 1 onward to randomize positions
  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(atomIndexVector.begin() + 1, atomIndexVector.end(), gen);

  tiledSupercell.UpdateAtomVector(atomIndexVector);

  auto vacancyId = tiledSupercell.GetVacancySiteId();

  cout << vacancyId.latticeId << " : " << vacancyId.smallConfigId << endl;

  tiledSupercell.WriteAtomicIndicesToFile("random_atomicIndices_primSC10_CS30_Mo50Ta50X.bin.gz");
}
  */