#include "Config.h"
#include "LVFEPredictor.h"
#include "SymmetricCEPredictor.h"

#include "ComputeVFE.h"


int main()
{

  size_t smallConfigSize = 10;
  double latticeParam = 3.2;
  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");
  size_t cubeSize = 35;
  Cube cubeObj(cubeSize);

  // 350x350x350 

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

  // Build the atomIndices vector from equilibriated configurations
  size_t numSites = tiledSupercell.GetTotalNumOfSites();


  std::vector<uint64_t> atomIndexVector(numSites);

  size_t finalIndex = 65700000;
  size_t stepsWidth = 100000;

  // 7x7x7 cube will be required to make a tiled supercell with 50x50x50 supercell

  string pathToConfigs = "//media/sf_Phd/MoTa/CMC/cmc1700K/configs"; 

  size_t idx = 0;

  for (int i = 0; i < 343; i++)
  {
    size_t currentIndex = finalIndex - stepsWidth*i;

    cout << currentIndex << endl;

    auto cfg = Config::ReadCfg(pathToConfigs + "/" + to_string(currentIndex) + ".cfg.gz");


    for (const auto ele : cfg.GetAtomVector())
    {
      atomIndexVector[idx] = uint64_t(ele.GetAtomicIndex());
      idx++;
    }
  }




  cout << "Number of sites: " << numSites << endl;
  cout << 250000*343 << endl;
  cout << idx << endl;

  string outPath = "//media/sf_Phd/MoTa/CMC/ss350x350x350/start_atomicIndices_primSC10_CS35_Mo50Ta50.bin.gz";

  tiledSupercell.UpdateAtomVector(atomIndexVector);

  tiledSupercell.WriteAtomicIndicesToFile(outPath);
  
}

/*
//// AVERAGE VACANCY CONCENTRATION
int main()
{
  string pathToCoefficients = "/home/pravendra3/Documents/TLMC/bin/coefficientFile_MoTa_V3.2.json";
  double latticeParam = 3.2;

  size_t supercellSize = 10;

  vector<double> cutoffs = {3, 4, 5};

  ClusterExpansionParameters ceParams(pathToCoefficients);

  auto supercellConfig = Config::GenerateAlloySupercell(
      supercellSize,
      latticeParam,
      "BCC",
      {"Mo", "Ta"},
      {50.0, 50.0},
      1);

  double maxClusterCutoff = ceParams.GetMaxClusterCutoff();
  supercellConfig.UpdateNeighborList({maxClusterCutoff});

  // Generate a small config for declaring symCE
  const size_t primSize = 2;

  auto primConfig = Config::GenerateSupercell(
      primSize,
      latticeParam,
      "Mo", // Does not matter as the lattice param and structure type is important
      "BCC");

  // Declare Symmetric CE
  SymmetricCEPredictor symCEEnergyPredictor(
      ceParams,
      supercellConfig,
      primConfig);

  // Again update the neighbor list
  supercellConfig.UpdateNeighborList(cutoffs);

  // Declare LVFE Predictor
  LVFEPredictor lvfePredictor(
      ceParams,
      supercellConfig);

  // Open the output file
  std::ofstream outFile("vfe_output.txt");

  // Get elements from chemical potential map for header
  auto chemicalPotMap = ceParams.GetChemicalPotentialsMap();
  std::vector<std::string> elementNames;
  for (const auto &entry : chemicalPotMap)
  {
    elementNames.push_back(entry.first);
  }

  // Write header
  outFile << "idx";
  for (const auto &eleName : elementNames)
  {
    outFile << "\t" << eleName;
  }
  outFile << "\n";
  outFile.flush();

  // Loop over iterations
  size_t numOfIterations = 100;
  for (size_t i = 0; i < numOfIterations; i++)
  {
    cout << "Iteration " << i << endl;
    auto currentConfig = Config::GenerateAlloySupercell(
        supercellSize,
        latticeParam,
        "BCC",
        {"Mo", "Ta"},
        {50.0, 50.0},
        i);

    auto vfeMap = ComputeVFE(
        currentConfig,
        chemicalPotMap,
        symCEEnergyPredictor,
        lvfePredictor);

    // Write index
    outFile << i;

    // Write VFE for each element in the same order as header
    for (const auto &eleName : elementNames)
    {
      outFile << "\t" << std::setprecision(8) << vfeMap.at(Element(eleName));
    }
    outFile << "\n";
    outFile.flush();
  }

  outFile.close();
}

*/
//// MAKING RANDOM CONFIG FOR TILED LATTICE MONTE CARLO SIMULATION
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
  size_t cubeSize = 5;
  Cube cubeObj(cubeSize);

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

  size_t numSites = tiledSupercell.GetTotalNumOfSites();

  // make a atom vector


  // Prepare vector
  std::vector<uint64_t> atomIndexVector(numSites);

  uint64_t moIndex = Element("Mo").GetAtomicIndex();
  uint64_t taIndex = Element("Ta").GetAtomicIndex();

  // Fill remaining sites
  size_t half = (numSites - 1) / 2; // minus first site
  for (size_t i = 0; i <= half; i++)
    atomIndexVector[i] = moIndex;
  for (size_t i = half + 1; i < numSites; i++)
    atomIndexVector[i] = taIndex;

  // Shuffle from index 1 onward to randomize positions
  std::random_device rd;
  std::mt19937 gen(rd());
  std::shuffle(atomIndexVector.begin() + 1, atomIndexVector.end(), gen);

  tiledSupercell.UpdateAtomVector(atomIndexVector);

  string filename = "random_atomicIndices_primSC"  + to_string(smallConfigSize)
  + "_CS" + to_string(cubeSize) + "_Mo50Ta50.bin.gz";
  tiledSupercell.WriteAtomicIndicesToFile(filename);
}
*/