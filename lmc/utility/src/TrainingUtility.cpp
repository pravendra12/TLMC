#include "TrainingUtility.h"

VectorXd GetKRAEncoding(
    const Config &config,
    const pair<size_t, size_t> &latticeIdPair,
    BasisSet &atomicBasis,
    const size_t &maxBondOrder,
    const size_t &maxClusterSize,
    const unordered_map<size_t, Eigen::RowVector3d> &referenceLatticeIdHashmap,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &encodedOrbits,
    const vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> &symmetryOperations)
{

  auto sortedLatticeIds = GetCanonicalSortedSitesForPair(
      config,
      latticeIdPair,
      maxBondOrder,
      referenceLatticeIdHashmap,
      symmetryOperations);

  VectorXd correlationVector = GetCorrelationVector(
      config,
      atomicBasis,
      sortedLatticeIds,
      encodedOrbits); 

  return correlationVector;
}

// Returns local cluster expansion around a lattice site
// To define local environment around a site which is
// Used for defining the local vacancy formation energy
VectorXd GetLocalSiteClusterVector(
    const Config &config,
    const size_t &latticeId,
    BasisSet &atomicBasis,
    const size_t maxBondOrder,
    const size_t maxClusterSize,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &encodedOrbits)
{
  // Canonically Sorted Lattice Ids around a site
  auto sortedIdsAroundSite = GetCanonicalSortedSitesForSite(
      config,
      latticeId,
      maxBondOrder);

  // Defines correlation or cluster vector around a lattice site to define
  // local cluster expansion
  // Uses ising type variable
  // For binary is valid, need to work on this for ternary and quaternary
  // A = -1 B = 1
  VectorXd correlationVector = GetCorrelationVector(
      config,
      atomicBasis,
      sortedIdsAroundSite,
      encodedOrbits); // so as to avoid a redundant constant term

  return correlationVector;
}

json ExtractTrainingDataForNEB(
    const string &pathToNEBOutput)
{
  const vector<double> cutoffs = {3, 4, 5};

  auto cfg = Config::GenerateSupercell(4, 3.4, "Mo", "BCC");
  cfg.UpdateNeighborList(cutoffs);

  const size_t maxClusterSize = 3;
  const size_t maxBondOrder = 3;
  const string basisType = "Chebyshev";

  set<Element> elementSet = {Element("Mo"), Element("Ta")};
  BasisSet atomicBasis(elementSet, basisType);

  Vector3d referenceJumpDirection(1, 1, 1);
  auto referenceLatticeIdHashmap = GetCenteredNeighborsAlongJumpDirection(
      cfg, maxBondOrder, referenceJumpDirection);

  auto centralLatticeId = cfg.GetCentralAtomLatticeId();
  auto centralPos = cfg.GetRelativePositionOfLattice(centralLatticeId);

  // Target direction <111>
  Vector3d target(1, 1, 1);
  target.normalize();

  size_t chosenNeighbor = -1; // placeholder
  double bestDot = -2.0;      // cosine similarity score

  // Search over neighbors
  for (int order = 1; order <= maxBondOrder; ++order)
  {
    auto neighbors = cfg.GetNeighborLatticeIdVectorOfLattice(centralLatticeId, order);
    for (auto neighborId : neighbors)
    {
      auto neighborPos = cfg.GetRelativePositionOfLattice(neighborId);
      Vector3d relVec = neighborPos - centralPos;
      relVec.normalize();

      double dot = relVec.dot(target);

      if (dot > bestDot) // best alignment with (1,1,1)
      {
        bestDot = dot;
        chosenNeighbor = neighborId;
      }
    }
  }
  cout << "Chosen neighbor = " << chosenNeighbor
       << " aligned with <111>, dot = " << bestDot << endl;

  // Now use this neighbor in your symmetry operation setup:
  auto nnLatticeIdSet = cfg.GetNeighboringLatticeIdSetOfPair(
      {centralLatticeId, chosenNeighbor},
      maxBondOrder);

  auto symOpsBCC = GetSymmetryOperations(cfg);

  auto encodedOrbitsForPair = GetLocalEncodedOrbitsForPair(cfg, maxBondOrder, maxClusterSize,
                                                           referenceJumpDirection, true);

  auto encodedOrbitsForSite = GetLocalEncodedOrbitsForSite(cfg, maxBondOrder, maxClusterSize, true);

  elementSet.insert(Element("X"));
  ConfigEncoding configEncoder(cfg, elementSet, maxBondOrder, maxClusterSize);

  json outputJson;

  // Iterate over each NEB path folder
  for (const auto &entry : fs::directory_iterator(pathToNEBOutput))
  {
    if (!entry.is_directory())
      continue;

    string folderPath = entry.path().string();
    string structuresPath = folderPath + "/structures/unrelaxed/";

    string poscarInitial = structuresPath + "POSCAR_unrelaxed_initial";
    string poscarFinal = structuresPath + "POSCAR_unrelaxed_final";

    if (!fs::exists(poscarInitial) || !fs::exists(poscarFinal))
    {
      cerr << "POSCAR files not found in " << structuresPath << endl;
      continue;
    }

    auto cfgInitial = Config::ReadPoscar(poscarInitial);
    auto cfgFinal = Config::ReadPoscar(poscarFinal);

    cfgInitial.UpdateNeighborList(cutoffs);
    cfgFinal.UpdateNeighborList(cutoffs);

    // Optional: neighbor sanity check
    auto check_neighbors = [](const Config &cfg, const string &label)
    {
      bool ok = true;
      if (cfg.GetNeighborLatticeIdVectorOfLattice(0, 1).size() != 8)
      {
        cerr << label << ": 1st neighbors != 8\n";
        ok = false;
      }
      if (cfg.GetNeighborLatticeIdVectorOfLattice(0, 2).size() != 6)
      {
        cerr << label << ": 2nd neighbors != 6\n";
        ok = false;
      }
      if (cfg.GetNeighborLatticeIdVectorOfLattice(0, 3).size() != 12)
      {
        cerr << label << ": 3rd neighbors != 12\n";
        ok = false;
      }
      if (!ok)
        throw runtime_error(label + ": neighbor check failed!");
    };
    check_neighbors(cfgInitial, "cfgInitial");
    check_neighbors(cfgFinal, "cfgFinal");

    // Vacancy pair
    auto latticeIdJumpPair = make_pair(cfgInitial.GetVacancyLatticeId(),
                                       cfgFinal.GetVacancyLatticeId());

    // Encode initial and final configs
    VectorXd ceEncodingInitial = configEncoder.GetEncodeVector(cfgInitial);
    VectorXd localClusterEncodingInitial = GetLocalSiteClusterVector(
        cfgInitial, latticeIdJumpPair.first, atomicBasis,
        maxBondOrder, maxClusterSize, encodedOrbitsForSite);

    VectorXd ceEncodingFinal = configEncoder.GetEncodeVector(cfgFinal);
    VectorXd localClusterEncodingFinal = GetLocalSiteClusterVector(
        cfgFinal, latticeIdJumpPair.second, atomicBasis,
        maxBondOrder, maxClusterSize, encodedOrbitsForSite);

    // KRA encoding
    VectorXd kraEncoding = GetKRAEncoding(
        cfgInitial, latticeIdJumpPair, atomicBasis, maxBondOrder, maxClusterSize,
        referenceLatticeIdHashmap, encodedOrbitsForPair, symOpsBCC);

    // Store in JSON
    json entryJson;
    string folderName = entry.path().filename().string(); // <-- base folder name
    entryJson["folder"] = folderName;
    entryJson["jumpPair"] = {latticeIdJumpPair.first, latticeIdJumpPair.second};
    entryJson["migratingElement"] = cfgInitial.GetElementOfLattice(latticeIdJumpPair.second).GetElementString();
    entryJson["ceEncodingInitial"] = vector<double>(ceEncodingInitial.data(),
                                                    ceEncodingInitial.data() + ceEncodingInitial.size());
    entryJson["localClusterEncodingInitial"] = vector<double>(localClusterEncodingInitial.data(),
                                                              localClusterEncodingInitial.data() + localClusterEncodingInitial.size());
    entryJson["ceEncodingFinal"] = vector<double>(ceEncodingFinal.data(),
                                                  ceEncodingFinal.data() + ceEncodingFinal.size());
    entryJson["localClusterEncodingFinal"] = vector<double>(localClusterEncodingFinal.data(),
                                                            localClusterEncodingFinal.data() + localClusterEncodingFinal.size());
    entryJson["kraEncoding"] = vector<double>(kraEncoding.data(),
                                              kraEncoding.data() + kraEncoding.size());

    outputJson[entry.path().filename().string()] = entryJson;

    cout << "Done for " << folderName << endl;
  }

  return outputJson;
}

json ExtractLocalCEData(const string &pathToLocalCEOutput)
{
  const vector<double> cutoffs = {3, 4, 5};

  auto cfg = Config::GenerateSupercell(4, 3.2, "Mo", "BCC");
  cfg.UpdateNeighborList(cutoffs);

  const size_t maxClusterSize = 3;
  const size_t maxBondOrder = 3;
  const string basisType = "Chebyshev";

  set<Element> elementSet = {Element("Mo"), Element("Ta")};
  BasisSet atomicBasis(elementSet, basisType);

  auto centralLatticeId = cfg.GetCentralAtomLatticeId();
  auto encodedOrbitsForSite = GetLocalEncodedOrbitsForSite(cfg, maxBondOrder, maxClusterSize, true);

  elementSet.insert(Element("X"));
  ConfigEncoding configEncoder(cfg, elementSet, maxBondOrder, maxClusterSize);

  json outputJson;

  for (const auto &entry : fs::directory_iterator(pathToLocalCEOutput))
  {
    if (!entry.is_directory())
      continue;

    string folderPath = entry.path().string();
    string folderName = entry.path().filename().string(); // base folder name

    string moPath = folderPath + "/Mo/structure/POSCAR_unrelaxed_Mo";
    if (!fs::exists(moPath))
    {
      cerr << "POSCAR_unrelaxed_Mo not found in " << folderPath << endl;
      continue;
    }

    auto cfgMo = Config::ReadPoscar(moPath);
    cfgMo.UpdateNeighborList(cutoffs);

    // Optional: sanity check
    auto check_neighbors = [](const Config &cfg, const string &label)
    {
      if (cfg.GetNeighborLatticeIdVectorOfLattice(0, 1).size() != 8 ||
          cfg.GetNeighborLatticeIdVectorOfLattice(0, 2).size() != 6 ||
          cfg.GetNeighborLatticeIdVectorOfLattice(0, 3).size() != 12)
      {
        throw runtime_error(label + ": neighbor check failed!");
      }
    };
    check_neighbors(cfgMo, folderName);

    // Prepare JSON entry
    json entryJson;

    // Mo
    entryJson["Mo"]["ceEncoding"] = configEncoder.GetEncodeVector(cfgMo);
    entryJson["Mo"]["localClusterEncoding"] = GetLocalSiteClusterVector(
        cfgMo, centralLatticeId, atomicBasis,
        maxBondOrder, maxClusterSize, encodedOrbitsForSite);

    // Ta
    cfgMo.SetElementOfLattice(centralLatticeId, Element("Ta"));
    entryJson["Ta"]["ceEncoding"] = configEncoder.GetEncodeVector(cfgMo);
    entryJson["Ta"]["localClusterEncoding"] = GetLocalSiteClusterVector(
        cfgMo, centralLatticeId, atomicBasis,
        maxBondOrder, maxClusterSize, encodedOrbitsForSite);

    // X (vacancy)
    cfgMo.SetElementOfLattice(centralLatticeId, Element("X"));
    entryJson["X"]["ceEncoding"] = configEncoder.GetEncodeVector(cfgMo);
    entryJson["X"]["localClusterEncoding"] = GetLocalSiteClusterVector(
        cfgMo, centralLatticeId, atomicBasis,
        maxBondOrder, maxClusterSize, encodedOrbitsForSite);

    // Add to main JSON
    outputJson[folderName] = entryJson;

    cout << "Done for " << folderName << endl;
  }

  return outputJson;
}

void IterateDirectoryToExtractData()
{
  string pathToNebOutput0 = "/home/pravendra3/Documents/nebOutput/nebOutput";
  string pathToNebOutput1 = "/home/pravendra3/Documents/nebOutput/activeLearning1";
  string pathToNebOutput2 = "/home/pravendra3/Documents/nebOutput/activeLearning2";
  string pathToNebOutput3 = "/home/pravendra3/Documents/nebOutput/activeLearning3";
  string pathToNebOutput4 = "/home/pravendra3/Documents/nebOutput/activeLearning4";

  // Collect NEB outputs
  json nebOutput;
  nebOutput.push_back(ExtractTrainingDataForNEB(pathToNebOutput0));
  nebOutput.push_back(ExtractTrainingDataForNEB(pathToNebOutput1));
  nebOutput.push_back(ExtractTrainingDataForNEB(pathToNebOutput2));
  nebOutput.push_back(ExtractTrainingDataForNEB(pathToNebOutput3));
  nebOutput.push_back(ExtractTrainingDataForNEB(pathToNebOutput4));

  // Write combined NEB output to file
  {
    std::ofstream ofs("/home/pravendra3/Documents/nebOutput/encodingForNEBData.json");
    ofs << nebOutput.dump(2); // pretty print with indent=2
  }

  // Local CE
  string pathToLocalCEOutput = "/home/pravendra3/Documents/nebOutput/relaxationOutput";
  json lceOutput = ExtractLocalCEData(pathToLocalCEOutput);

  // Write Local CE output to file
  {
    std::ofstream ofs("/home/pravendra3/Documents/nebOutput/encodingForLocalCEData.json");
    ofs << lceOutput.dump(2);
  }
}