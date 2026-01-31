
#include "Home.h"

/*
int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    std::cout << "No input parameter filename." << std::endl;
    return 1;
  }
  api::Parameter parameter(argc, argv);
  api::Print(parameter);
  api::Run(parameter);
}
*/

#include <filesystem>
#include <string>
#include <iostream>

namespace fs = std::filesystem;

int main()
{

  // CubeSize, supercellSize
  vector<pair<size_t, size_t>> supercellSizeVector = {
      {3, 10},
  };

  const string pathProfileSupercellSize = "/home/pravendra3/Documents/TLMC/bin/profilingSupercellSize";

  const double latticeParam = 3.2;

  for (const auto &cubeAndSCSize : supercellSizeVector)
  {
    const size_t smallConfigSize = cubeAndSCSize.second;
    const size_t cubeSize = cubeAndSCSize.first;

    auto smallCfg = Config::GenerateSupercell(
        smallConfigSize,
        latticeParam,
        "Mo",
        "BCC");
    Cube cubeObj(cubeSize);

    TiledSupercell tiledSupercell(
        smallCfg,
        cubeObj);

    size_t numSites = tiledSupercell.GetTotalNumOfSites();

    auto largeConfig = Config::GenerateSupercell(
        30,
        3.2,
        "Mo",
        "BCC");

    for (size_t i = 0; i < numSites / 2; i++)
    {
      largeConfig.SetElementOfLattice(i, Element("Ta"));
    }

    largeConfig.SetElementOfLattice(largeConfig.GetCentralAtomLatticeId(), Element("X"));

    string currentDir = pathProfileSupercellSize + "/primSC" + to_string(smallConfigSize) + "CS" + to_string(cubeSize);

    fs::path dirPath(currentDir);

    // Create the directory if it doesn't exist
    if (!fs::exists(dirPath))
    {
      fs::create_directories(dirPath);
      std::cout << "Created directory: " << dirPath << "\n";
    }

    // Generate filenames
    std::string atomicIndexFilename = "random_atomicIndices_primSC" + std::to_string(smallConfigSize) +
                                      "_CS" + std::to_string(cubeSize) + "_Mo50Ta50";
    std::string configFilename = atomicIndexFilename;

    // Full paths
    fs::path atomicIndexFilePath = dirPath / (atomicIndexFilename + ".bin.gz");
    fs::path configFilePath = dirPath / (configFilename + ".cfg.gz");
    Config::WriteConfig(configFilePath.string(), largeConfig);

    auto atomVector = largeConfig.GetAtomVector();

    // make a atom vector

    // Prepare vector
    std::vector<uint64_t> atomIndexVector(numSites);

    for (size_t i = 0; i < numSites; i++)
    {
      atomIndexVector[i] = atomVector[i].GetAtomicIndex();
    }

    tiledSupercell.UpdateAtomVector(atomIndexVector);

    // Write atomic indices
    tiledSupercell.WriteAtomicIndicesToFile(atomicIndexFilePath.string());

    std::cout << "Files written to: " << dirPath << "\n";
  }
}
