#include "ConvertAtomVectorsToConfigs.h"

void ConvertAtomVectorsToConfigs(const string &pathToTLMCRun)
{
  cout << "[INFO] Starting conversion in: " << pathToTLMCRun << endl;

  string paramFileCmc = pathToTLMCRun + "/cmc_param.txt";
  string paramFileKmc = pathToTLMCRun + "/kmc_param.txt";

  string paramFile;

  if (fs::exists(paramFileCmc))
  {
    paramFile = paramFileCmc;
    cout << "[INFO] Using parameter file: cmc_param.txt" << endl;
  }
  else if (fs::exists(paramFileKmc))
  {
    paramFile = paramFileKmc;
    cout << "[INFO] Using parameter file: kmc_param.txt" << endl;
  }
  else
  {
    cerr << "Error in `ConvertAtomVectorsToConfigs`: Expected either 'cmc_param.txt' or 'kmc_param.txt' in " << pathToTLMCRun << endl;
    return;
  }

  api::Parameter parameter(paramFile);
  cout << "[INFO] Loaded parameters successfully." << endl;

  // Generate small configuration and cube
  auto smallCfg = Config::GenerateSupercell(
      parameter.supercell_size_,
      parameter.lattice_param_,
      "X", // placeholder element
      parameter.structure_type_);
  cout << "[INFO] Generated small supercell (placeholder X)." << endl;

  Cube cubeObj(parameter.cube_size_);
  TiledSupercell tiledSupercell(smallCfg, cubeObj);
  cout << "[INFO] Created TiledSupercell object." << endl;

  // Ensure output directory exists
  const string outputPath = pathToTLMCRun + "/configs";
  fs::create_directories(outputPath);
  cout << "[INFO] Ensured output directory exists: " << outputPath << endl;

  // Loop through all files in the folder
  for (const auto &entry : fs::directory_iterator(pathToTLMCRun))
  {
    auto filePath = entry.path();

    // Process only .bin.gz files
    if (filePath.extension() != ".gz" || filePath.stem().extension() != ".bin")
      continue;

    string baseName = filePath.stem().stem().string(); // "0" from "0.bin.gz"
    cout << "[INFO] Processing file: " << filePath.filename() << endl;

    // Read atom indices
    auto atomIndicesVector = TiledSupercell::ReadAtomicIndicesFromFile(filePath.string());

    // Convert to Element objects
    vector<Element> atomVector;
    atomVector.reserve(atomIndicesVector.size());
    for (auto idx : atomIndicesVector)
      atomVector.emplace_back(Element(idx));

    // Build large config
    auto largeCfgFromAtomVector = TiledSupercell::MakeSupercellFromAtomInfo(
        smallCfg, cubeObj, atomVector);

    // Write output
    string outFile = outputPath + "/" + baseName + ".cfg.gz";
    Config::WriteConfig(outFile, largeCfgFromAtomVector);
  }

  cout << "[INFO] Conversion completed successfully." << endl;
}
