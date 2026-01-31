/*! \file  main.cpp
 *  \brief File for the main function.
 */

#include "Home.h"
#include <filesystem>
#include <string>
#include <iostream>
#include "SaveCubeAsConfig.h"

namespace fs = std::filesystem;

vector<size_t> GetSelectedCubeIds(const string &filepath)
{
  // read the stored cubeIds
  std::ifstream fin(filepath);
  if (!fin)
  {
    std::cerr << "Failed to open file\n";
  }

  // skip first two lines
  // num of cubes
  // lattice info
  // do not do it now

  std::vector<size_t> selectedCubeIds;
  std::string line;

  while (std::getline(fin, line))
  {
    std::stringstream ss(line);
    std::string token;

    while (std::getline(ss, token, ',')) // handles comma separation
    {
      try
      {
        size_t val = static_cast<size_t>(std::stoull(token));
        // cout << val << endl;
        selectedCubeIds.push_back(val);
      }
      catch (...)
      {
        std::cerr << "❌ Failed token: [" << token << "]\n";
      }
    }
  }
  return selectedCubeIds;
}

Config GetSlicedConfig(
    const TiledSupercell &tiledSupercell,
    const vector<size_t> &selectedCubeIds,
    const size_t numOfCubeX = 10,
    const size_t numOfCubeY = 10)
{

  size_t numOfCubes = (selectedCubeIds.size() - 1); // first entry is the total number of cubes
  size_t numOfCubeZ = size_t(numOfCubes / (numOfCubeX * numOfCubeY));
  size_t numOfSitesPerSmallConfig = tiledSupercell.GetNumOfSitesPerSmallConfig();

  Matrix3d basisOfSmallConfig = tiledSupercell.GetSmallConfig().GetBasis();

  Matrix3d basisOfSlicedConfig;
  basisOfSlicedConfig.col(0) = Vector3d({numOfCubeX * basisOfSmallConfig(0, 0), 0, 0});
  basisOfSlicedConfig.col(1) = Vector3d({0, numOfCubeY * basisOfSmallConfig(1, 1), 0});
  basisOfSlicedConfig.col(2) = Vector3d({0, 0, numOfCubeZ * basisOfSmallConfig(2, 2)});

  // For cubes or smallConfigs lying inside the sliced region
  Matrix3Xi reducedRelativeIndexMatrix(3, numOfCubes);

  for (size_t i = 1; i < selectedCubeIds.size(); i++)
  {
    size_t cubeId = selectedCubeIds[i];
    reducedRelativeIndexMatrix.col(i - 1) = tiledSupercell.GetCube().GetRelativePosition(cubeId);
  }

  Vector3i minVector = reducedRelativeIndexMatrix.rowwise().minCoeff();
  // Shift all columns so that the minimum becomes zero
  reducedRelativeIndexMatrix.colwise() -= minVector;

  vector<Element> originalAtomVector = tiledSupercell.GetAtomVector();
  const Matrix3Xd relativePositionMatrixSmallCfg = tiledSupercell.GetSmallConfig().GetRelativePositionMatrix();

  // Sliced Config
  size_t numOfAtomsSlicedConfig = numOfSitesPerSmallConfig * numOfCubes;

  Matrix3Xd relativePosMatrixSlicedConfig(3, numOfAtomsSlicedConfig);

  vector<Element> atomVectorSlicedConfig;
  atomVectorSlicedConfig.reserve(numOfAtomsSlicedConfig);

  size_t colIdx = 0;
  for (size_t cubeId = 0; cubeId < numOfCubes; ++cubeId)
  {
    // 1) Get offset for this small cube, in cube units
    Vector3d offSet = reducedRelativeIndexMatrix.col(cubeId).cast<double>();

    // 2) Add offset in cube space
    Matrix3Xd transformedRelativePositions =
        relativePositionMatrixSmallCfg.colwise() + offSet;

    // 3) Convert from "cube units" to fractional coords of the new supercell
    transformedRelativePositions.row(0) /= static_cast<double>(numOfCubeX);
    transformedRelativePositions.row(1) /= static_cast<double>(numOfCubeY);
    transformedRelativePositions.row(2) /= static_cast<double>(numOfCubeZ);

    // 4) Insert block directly into the sliced config matrix
    relativePosMatrixSlicedConfig.middleCols(
        colIdx,
        numOfSitesPerSmallConfig) = transformedRelativePositions;

    // 5) Append the corresponding elements for this cube
    size_t smallConfigIdx = selectedCubeIds[cubeId + 1]; // first entry is the total number of cubes

    size_t startIndex = smallConfigIdx * numOfSitesPerSmallConfig;
    size_t endIndex = startIndex + numOfSitesPerSmallConfig;

    atomVectorSlicedConfig.insert(atomVectorSlicedConfig.end(),
                                  originalAtomVector.begin() + startIndex,
                                  originalAtomVector.begin() + endIndex);

    colIdx += numOfSitesPerSmallConfig;
  }

  return Config(
      basisOfSlicedConfig,
      relativePosMatrixSlicedConfig,
      atomVectorSlicedConfig);
}

int main()
{
  string pathTiledSupercellOutput = "//media/sf_Phd/TiledLMC";

  string pathTLMCOutput = "//media/sf_Phd/MoTa/run_tlmc/ss350x350x350/ceMoTaX/data/raw/qSim/cc100um/1400K_to_300K";
  string pathSimulateTEM = "//media/sf_Phd/MoTa/run_tlmc/ss350x350x350/ceMoTaX/data/results/simulated_diffraction_pattern/withTDS/sliced_bcc/qSim/cc100um/1400K_to_300K";

  vector<string> pathToSelectedCubeIdFileVector = {
      "sliced_z_30nm_001_ZA_view_cube_Ids",
      "sliced_z_45nm_001_ZA_view_cube_Ids",
      "sliced_z_60nm_001_ZA_view_cube_Ids",
  };

  size_t cubeSize = 35;
  size_t smallConfigSize = 10;
  const double latticeParam = 3.2;

  size_t numOfCubeX = 10;
  size_t numOfCubeY = 10;

  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");

  const size_t numLatticePerConfig = smallCfg.GetNumLattices();

  Cube cubeObj(cubeSize);

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

//  for (const auto &entry : fs::directory_iterator(pathTLMCOutput))
//  {
//    if (!entry.is_directory())
//      continue;
//
//    string simFolder = entry.path().filename().string();
//
//    // if (!(simFolder.starts_with("cmc") || simFolder.starts_with("kmc")))
//    //   continue;
//
//    std::cout << "\n📁 Processing simulation folder: "
//              << simFolder << "\n";
//
//    string simFolderPath = entry.path().string();

    // -------------------------------------------------------
    // ITERATE OVER *ALL* .bin.gz FILES IN THIS FOLDER
    // -------------------------------------------------------
    for (const auto &binEntry : fs::directory_iterator(pathTLMCOutput))
    {
      string filename = binEntry.path().filename().string();

      if (!(filename.ends_with(".bin.gz")))
        continue;

      string binFilePath = binEntry.path().string();

      // folder name = just timestep number
      string binBaseName = filename;
      size_t pos = binBaseName.find(".bin.gz");
      if (pos != string::npos)
        binBaseName = binBaseName.substr(0, pos);

      std::cout << "   ✅ Found bin: " << binBaseName << "\n";

      // -------------------------------------------------------
      // LOAD ATOMIC CONFIG FROM BIN
      // -------------------------------------------------------
      auto atomicIndexVector =
          TiledSupercell::ReadAtomicIndicesFromFile(binFilePath);

      tiledSupercell.UpdateAtomVector(atomicIndexVector);

      // -------------------------------------------------------
      // ITERATE OVER ALL SLICE DEFINITIONS
      // -------------------------------------------------------
      for (const auto &selectedCubeIdFile : pathToSelectedCubeIdFileVector)
      {
        string pathToSelectedCubeIdFile =
            pathTiledSupercellOutput + "/" + selectedCubeIdFile;

        auto selectedCubeIds =
            GetSelectedCubeIds(pathToSelectedCubeIdFile);

        // strip suffix for folder name
        string sliceName = selectedCubeIdFile;
        size_t pos2 = sliceName.find("_cube_Ids");
        if (pos2 != string::npos)
          sliceName = sliceName.substr(0, pos2);

        // -------------------------------------------------------
        // ✅ FINAL OUTPUT DIRECTORY
        // SimulateTEM/cmc1000K/1000000000/sliced_z_30nm_001_ZA/
        // -------------------------------------------------------
        // string outFolder =
        //     pathSimulateTEM + "/" +
        //     simFolder + "/" +
        //     binBaseName + "/" +
        //     sliceName;


        string outFolder =
            pathSimulateTEM + "/" +
            binBaseName + "/" +
            sliceName;

        // string outFolder = pathSimulateTEM;

        fs::create_directories(outFolder);

        // -------------------------------------------------------
        // BUILD SLICED CONFIG
        // -------------------------------------------------------

        // -------------------------------------------------------
        // WRITE RAW CONFIG
        // -------------------------------------------------------
        size_t numOfCubes = (selectedCubeIds.size() - 1);
        size_t numAtoms = numLatticePerConfig * numOfCubes;

        vector<Element> atomVector;
        atomVector.reserve(numAtoms);

        Matrix3Xd relativePositionMatrix(3, numAtoms);

        size_t colIdx = 0;

        for (size_t i = 1; i < selectedCubeIds.size(); i++)
        {
          size_t cubeId = selectedCubeIds[i];

          for (size_t latticeId = 0;
               latticeId < smallCfg.GetNumLattices();
               latticeId++)
          {
            LatticeSiteMapping latticeSite;
            latticeSite.latticeId = latticeId;
            latticeSite.smallConfigId = cubeId;

            auto element =
                tiledSupercell.GetElementAtSite(latticeSite);

            atomVector.emplace_back(element);

            Vector3d relativePos =
                tiledSupercell.GetRelativePositionOfLatticeSiteMapping(latticeSite);

            relativePositionMatrix.col(colIdx) = relativePos;
            colIdx++;
          }
        }

        Config slicedConfig(
            tiledSupercell.GetSuperBasis(),
            relativePositionMatrix,
            atomVector);

        string rawOut =
            outFolder + "/sliced_raw.cfg.gz";

        Config::WriteConfig(rawOut, slicedConfig);

        // -------------------------------------------------------
        // WRITE TIGHT SUPERCELL CONFIG
        // -------------------------------------------------------
        Config slicedConfigShifted = GetSlicedConfig(
            tiledSupercell,
            selectedCubeIds);

        string tightOut =
            outFolder + "/sliced_tight.cfg.gz";

        Config::WriteConfig(tightOut, slicedConfigShifted);

        std::cout << "      ✅ Wrote slice: "
                  << sliceName << "\n";
      }
    }
  // }

  std::cout << "\n✅✅ ALL BIN FILES + SLICES PROCESSED ✅✅\n";
}
