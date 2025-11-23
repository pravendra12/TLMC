
#include "B2Cluster.h"
#include "B2ClusterTLMC.h"

#include <chrono>
using namespace std::chrono;
int main()
{

  size_t smallConfigSize = 10;
  double latticeParam = 3.2;
  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");

  smallCfg.UpdateNeighborList({3, 4});

  size_t cubeSize = 5;
  Cube cubeObj(cubeSize);

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

  auto atomVector = TiledSupercell::ReadAtomicIndicesFromFile("/home/pravendra3/Documents/TLMC/bin/testingB2ClusterTLMC/1000000.bin.gz");

  tiledSupercell.UpdateAtomVector(atomVector);

  tiledSupercell.UpdateNeighbourLists(2);

  // Time b2Cluster (TLMC version)
  auto t1 = high_resolution_clock::now();
  B2ClusterTLMC b2Cluster(tiledSupercell);
  auto t2 = high_resolution_clock::now();

  std::cout << "Time to create b2Cluster = "
            << duration_cast<milliseconds>(t2 - t1).count()
            << " ms\n";

  b2Cluster.WriteB2ClusterConfig(
      "/home/pravendra3/Documents/TLMC/bin/testingB2ClusterTLMC/tlmc_processed_1000000.xyz.gz");

  auto clusterMap = b2Cluster.GetB2Clusters();

  cout << clusterMap.size() << endl;

  auto largeConfig = tiledSupercell.MakeSupercell();
  largeConfig.UpdateNeighborList({3, 4});

  // Time b2ClusterOriginal creation
  auto t5 = high_resolution_clock::now();
  B2Cluster b2ClusterOriginal(largeConfig);
  auto t6 = high_resolution_clock::now();

  std::cout << "Time to create b2ClusterOriginal = "
            << duration_cast<milliseconds>(t6 - t5).count()
            << " ms\n";

  b2ClusterOriginal.WriteB2ClusterConfig(
      "/home/pravendra3/Documents/TLMC/bin/testingB2ClusterTLMC/cfg_processed_1000000.xyz.gz");

  auto clusters = b2ClusterOriginal.GetB2Clusters();
  cout << clusters.size() << endl;

}

/// SHORT RANGE ORDERING
/*
#include <filesystem>
#include "B2Cluster.h"
#include "B2ClusterTLMC.h"
#include "ShortRangeOrder.h"
#include "ShortRangeOrderTLMC.h"
#include <chrono>
using namespace std::chrono;
namespace fs = std::filesystem;

int main()
{
  size_t smallConfigSize = 10;
  double latticeParam = 3.2;

  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");

  const vector<double> cutoffs = {3.0};
  smallCfg.UpdateNeighborList(cutoffs);

  size_t cubeSize = 5;
  Cube cubeObj(cubeSize);

  TiledSupercell tiledSupercell(smallCfg, cubeObj);

  set<Element> elementSet = {Element("Mo"), Element("Ta"), Element("X")};

  string baseDir = "/media/sf_Phd/MoTa/validation/validation/MoTaX/cmc1000K_v4/";

  string configDir = baseDir + "/configs/";

  for (const auto &entry : fs::directory_iterator(configDir))
  {

    if (entry.path().extension() == ".gz" &&
        entry.path().string().find(".xyz") != string::npos)
    {
      string xyzFile = entry.path().string();

      // Corresponding atom vector file
      // Example: 13000000.xyz.gz  →  13000000.bin.gz
      string stem = entry.path().stem().stem().string(); // removes .xyz.gz
      string binFile = baseDir + "/" + stem + ".bin.gz";

      cout << "\n=============================\n";
      cout << "Processing configuration: " << stem << endl;
      cout << "XYZ file: " << xyzFile << endl;
      cout << "BIN file: " << binFile << endl;
      cout << "=============================\n";

      // --- Read large config ---
      auto largeConfig = Config::ReadXyz(xyzFile);
      largeConfig.UpdateNeighborList(cutoffs);

      // --- Compute SRO from large config ---
      ShortRangeOrder sro(largeConfig, elementSet);
      auto sroMap = sro.FindWarrenCowley(1);

      cout << "Warren-Cowley (Large Config):\n";
      for (const auto &entry : sroMap)
        cout << entry.first << " : " << entry.second << endl;

      // --- TLMC version ---
      auto atomVector = TiledSupercell::ReadAtomicIndicesFromFile(binFile);
      tiledSupercell.UpdateAtomVector(atomVector);
      tiledSupercell.UpdateNeighbourLists(cutoffs.size());

      ShortRangeOrderTLMC sroTLMC(tiledSupercell);
      auto eleProb = sroTLMC.ComputeWarrenCowley(1);

      cout << "Warren-Cowley (TLMC):\n";
      for (const auto &entry : eleProb)
        cout << entry.first << " : " << entry.second << endl;
    }
  }

  return 0;
}
*/

/// Testing Traverse

/*

#include "B2Cluster.h"
#include "B2ClusterTLMC.h"
#include "ShortRangeOrder.h"
#include "ShortRangeOrderTLMC.h"
#include <chrono>
#include "Traverse.h"
using namespace std::chrono;
int main()
{

  size_t smallConfigSize = 10;
  double latticeParam = 3.2;
  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");

  const vector<double> cutoffs = {3, 4};

  smallCfg.UpdateNeighborList(cutoffs);

  size_t cubeSize = 5;
  Cube cubeObj(cubeSize);

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

  tiledSupercell.UpdateNeighbourLists(cutoffs.size());

  ansys::Traverse iterator(
      0,
      1000000,
      "canonical_mc");

  SubLatticeOccupancy subLatticeOccupancy(
      tiledSupercell);

  set<Element> elementSet = {Element("Mo"), Element("Ta")};

  unordered_set<size_t> convertToConfigSet = {0, 1000000, 10000000};

  iterator.RunAnsys(
      tiledSupercell,
      subLatticeOccupancy,
      elementSet,
      convertToConfigSet);
}

*/