#include "TiledSupercell.h"
#include "Config.h"

TiledSupercell::TiledSupercell(
    const Config &smallConfig,
    const Cube &cubeObj) : smallConfig_(move(smallConfig)),
                           cubeObj_(move(cubeObj)),
                           numSitesPerSmallConfig_(
                               smallConfig.GetNumLattices()),
                           numSmallConfig_(
                               cubeObj.GetNumOfSites()),
                           totalNumOfSites_(
                               numSitesPerSmallConfig_ * numSmallConfig_)
{
  // Atom vector can be initialized and the neighbour lists can be updated after
  // the object has been declared based on the maxBondOrder // need to think about
  // it
  InitializeAtomVector();
  PrintTiledSupercell();
}

// these latticePair are bascially teh
// <LatticeId, smallConfigId> which will be used to
// get the atomId then the swap will happen
void TiledSupercell::LatticeJump(
    const pair<size_t, size_t> &latticeAndConfigIdPair1,
    const pair<size_t, size_t> &latticeAndConfigIdPair2)
{
  auto atomId1 = GetAtomIdFromLatticeAndConfigId(
      latticeAndConfigIdPair1);
  auto atomId2 = GetAtomIdFromLatticeAndConfigId(
      latticeAndConfigIdPair2);

  if (atomId1 >= atomVector_.size() || atomId2 >= atomVector_.size())
  {
    throw std::out_of_range(
        "Error in TiledSupercell::LatticeJump: atomId out of bounds "
        "(atomId1=" +
        std::to_string(atomId1) +
        ", atomId2=" + std::to_string(atomId2) +
        ", atomVector_.size()=" + std::to_string(atomVector_.size()) + ")");
  }

  std::swap(atomVector_[atomId1], atomVector_[atomId2]);
}

vector<pair<size_t, int>> TiledSupercell::GetNeighborLatticeIdVectorOfLattice(
    const size_t &latticeId,
    const size_t &distanceOrder) const
{
  if (distanceOrder == 0)
  {
    throw std::invalid_argument("Error in TiledSupercell::GetNeighborLatticeIdVectorOfLattice: distanceOrder must be >= 1");
  }
  if (distanceOrder > neighbourList_.size())
  {
    throw std::out_of_range(
        "Error in TiledSupercell::GetNeighborLatticeIdVectorOfLattice: distanceOrder (" +
        std::to_string(distanceOrder) +
        ") exceeds neighbourList_ size (" +
        std::to_string(neighbourList_.size()) +
        ")");
  }

  if (latticeId >= neighbourList_[distanceOrder - 1].size())
  {
    throw std::out_of_range(
        "Error in TiledSupercell::GetNeighborLatticeIdVectorOfLattice: latticeId (" +
        std::to_string(latticeId) +
        ") exceeds neighbourList_[distanceOrder-1] size (" +
        std::to_string(neighbourList_[distanceOrder - 1].size()) +
        ")");
  }

  // now safe to call at()
  return neighbourList_.at(distanceOrder - 1).at(latticeId);
}

void TiledSupercell::SetElementAtSite(
    const pair<size_t, size_t> &latticeAndConfigIdPair,
    const Element &element)
{
  size_t atomId = GetAtomIdFromLatticeAndConfigId(
      latticeAndConfigIdPair);

  if (atomId >= totalNumOfSites_)
  {
    throw std::out_of_range(
        "TiledSupercell::SetElementAtSite error: "
        "atomId " +
        std::to_string(atomId) +
        " (from latticeId " + std::to_string(latticeAndConfigIdPair.first) +
        ", smallConfigId " + std::to_string(latticeAndConfigIdPair.second) +
        ") exceeds total number of sites " + std::to_string(totalNumOfSites_));
  }
  atomVector_[atomId] = element;
}

Element TiledSupercell::GetElementAtSite(
    const pair<size_t, size_t> &latticeAndConfigIdPair)
{
  size_t atomId = GetAtomIdFromLatticeAndConfigId(
      latticeAndConfigIdPair);

  if (atomId >= totalNumOfSites_)
  {
    throw std::out_of_range(
        "TiledSupercell::GetElementAtSite error: "
        "atomId " +
        std::to_string(atomId) +
        " (from latticeId " + std::to_string(latticeAndConfigIdPair.first) +
        ", smallConfigId " + std::to_string(latticeAndConfigIdPair.second) +
        ") exceeds total number of sites " + std::to_string(totalNumOfSites_));
  }

  return atomVector_.at(atomId);
}

inline pair<size_t, size_t> TiledSupercell::GetLatticeAndConfigIdFromAtomId(
    const size_t &atomId) const
{
  if (atomId >= totalNumOfSites_)
  {
    throw std::out_of_range("Error in `TiledSupercell::GetLatticeAndConfigIdFromAtomId`: Atom index exceeds total number of atoms");
  }

  size_t smallConfigId = atomId / numSitesPerSmallConfig_;
  size_t latticeId = atomId % numSitesPerSmallConfig_;

  return {latticeId, smallConfigId};
}

inline size_t TiledSupercell::GetAtomIdFromLatticeAndConfigId(
    const pair<size_t, size_t> &latticeAndConfigIdPair) const
{
  size_t latticeId = latticeAndConfigIdPair.first;
  size_t smallConfigId = latticeAndConfigIdPair.second;

  size_t atomId = smallConfigId * numSitesPerSmallConfig_ + latticeId;

  if (atomId >= totalNumOfSites_)
  {
    throw std::out_of_range("Error in `TiledSupercell::GetAtomIdFromLatticeAndConfigId` : Computed atom index exceeds total number of atoms");
  }
  return atomId;
}

// I/O

Config TiledSupercell::MakeSupercell() const
{
  return MakeSupercellFromAtomInfo(
      smallConfig_,
      cubeObj_,
      atomVector_);
}

Config TiledSupercell::MakeSupercellFromAtomInfo(
    const Config &smallConfig,
    const Cube &cubeObj,
    const vector<Element> &atomVector)
{
  const size_t cubeSize = cubeObj.GetSizeOfCube();                      // returns a of cube
  const size_t numOfSitesPerSmallConfig = smallConfig.GetNumLattices(); // number of sites in each small config
  const size_t numOfSmallConfig = cubeObj.GetNumOfSites();
  const size_t totalNumOfSites = numOfSmallConfig * numOfSitesPerSmallConfig;

  Matrix3d superBasis = smallConfig.GetBasis() * static_cast<double>(cubeSize);

  Matrix3Xd relativePositionMatrixSmallCfg = smallConfig.GetRelativePositionMatrix(); // 3 x numAtomsSmallCfg
  Matrix3Xd relativePositionMatrixLargeCfg(3, totalNumOfSites);

  std::vector<Element> atomVectorLargeCfg;
  atomVectorLargeCfg.reserve(totalNumOfSites);

  size_t colIdx = 0;

  for (size_t smallConfigIdx = 0; smallConfigIdx < numOfSmallConfig; ++smallConfigIdx)
  {
    // 1) Get offset for this small cube
    Vector3d offSet = cubeObj.GetRelativePosition(smallConfigIdx).cast<double>();

    // 2) Add offset to *all* atoms of small config at once, scale by cubeSize
    Matrix3Xd transformedRelativePositions =
        (relativePositionMatrixSmallCfg.colwise() + offSet) / static_cast<double>(cubeSize);

    // 3) Insert block directly into the large matrix
    relativePositionMatrixLargeCfg.middleCols(
        colIdx,
        numOfSitesPerSmallConfig) = transformedRelativePositions;

    // 4) Append the corresponding elements for this cube
    size_t startIndex = smallConfigIdx * numOfSitesPerSmallConfig;
    size_t endIndex = startIndex + numOfSitesPerSmallConfig;

    atomVectorLargeCfg.insert(atomVectorLargeCfg.end(),
                              atomVector.begin() + startIndex,
                              atomVector.begin() + endIndex);

    colIdx += numOfSitesPerSmallConfig;
  }

  return Config(
      superBasis,
      relativePositionMatrixLargeCfg,
      atomVectorLargeCfg);
}

void TiledSupercell::WriteAtomVectorInfoToFile(
    const string &filename) const
{
  ofstream outFile(filename);
  if (!outFile.is_open())
  {
    throw std::runtime_error("Error in `TiledSupercell::WriteAtomVectorInfoToFile` : Could not open file: " + filename);
  }
  // optional header
  outFile << "IndexInAtomVector\tSmallCfgIdx\tLatticeId\tElementString\n";

  for (size_t i = 0; i < atomVector_.size(); ++i)
  {
    auto latticeAndConfigIdPair = GetLatticeAndConfigIdFromAtomId(i);

    size_t latticeId = latticeAndConfigIdPair.first;
    size_t smallCfgIdx = latticeAndConfigIdPair.second;

    outFile << i << "\t"
            << smallCfgIdx << "\t"
            << latticeId << "\t"
            << atomVector_[i].GetElementString() << "\n";
  }
}

void TiledSupercell::WriteAtomVectorInfoToBinary(
    const string &filename) const
{
  ofstream outFile(filename, ios::binary);
  if (!outFile)
    throw runtime_error("Error in `TiledSupercell::WriteAtomVectorInfoToBinary`: Could not open file");

  for (size_t i = 0; i < atomVector_.size(); ++i)
  {
    auto latticeAndConfigIdPair = GetLatticeAndConfigIdFromAtomId(i);

    size_t latticeId = latticeAndConfigIdPair.first;
    size_t smallCfgIdx = latticeAndConfigIdPair.second;

    size_t elementId = atomVector_[i].GetAtomicIndex(); // Atomic Number

    outFile.write(reinterpret_cast<const char *>(&smallCfgIdx), sizeof(size_t));
    outFile.write(reinterpret_cast<const char *>(&latticeId), sizeof(size_t));
    outFile.write(reinterpret_cast<const char *>(&elementId), sizeof(size_t));
  }
}

vector<Element> TiledSupercell::ReadAtomVectorInfoFromFile(
    const string &filename)
{
  ifstream inFile(filename);
  if (!inFile.is_open())
  {
    throw std::runtime_error("Error in `ReadAtomVectorInfoFromFile`: Could not open file: " + filename);
  }

  std::string headerLine;
  std::getline(inFile, headerLine); // skip header

  vector<Element> atomVector;

  size_t index, smallCfgIdx, latticeId;
  std::string elementStr;

  while (inFile >> index >> smallCfgIdx >> latticeId >> elementStr)
  {
    atomVector.push_back(Element(elementStr));
  }

  return atomVector;
}

vector<Element> TiledSupercell::ReadAtomVectorInfoFromBinary(
    const std::string &filename)
{
  std::ifstream inFile(filename, std::ios::binary);
  if (!inFile.is_open())
  {
    throw std::runtime_error("Error in `ReadAtomVectorInfoFromBinary`: Could not open file: " + filename);
  }

  vector<Element> atomVector;

  while (true)
  {
    size_t smallCfgIdx, latticeId, atomicNumber;

    // Read smallCfgIdx
    inFile.read(reinterpret_cast<char *>(&smallCfgIdx), sizeof(size_t));
    if (inFile.eof())
      break;

    // Read latticeId
    inFile.read(reinterpret_cast<char *>(&latticeId), sizeof(size_t));
    if (inFile.eof())
      break;

    // Read elementId
    inFile.read(reinterpret_cast<char *>(&atomicNumber), sizeof(size_t));
    if (inFile.eof())
      break;

    atomVector.emplace_back(Element(atomicNumber));
  }

  return atomVector;
}

// Helpers

// This function checks whether relative to a reference lattice Id the
// its neibhouringLatticeId whether it lies in same small config or it was
// neighbour due to PBC
bool TiledSupercell::LatticeInSameConfig(
    const size_t &referenceLatticeId,
    const size_t &neighbouringLatticeId) const
{
  const auto &basis = smallConfig_.GetBasis();

  // Distance with PBC
  auto relVecPBC = smallConfig_.GetRelativeDistanceVectorLattice(
      referenceLatticeId, neighbouringLatticeId);
  double distancePBC = (basis * relVecPBC).norm();

  // Distance without PBC
  auto diff = smallConfig_.GetRelativePositionOfLattice(referenceLatticeId) -
              smallConfig_.GetRelativePositionOfLattice(neighbouringLatticeId);
  double distanceNoPBC = (basis * diff).norm();

  // Compare with tolerance to handle floating-point issues
  return fabs(distancePBC - distanceNoPBC) < constants::kEpsilon;
}

// Finds which neighboring smallConfig (inside the cube) contains the neighbor site
// that is closest to the given reference site, considering periodic boundary conditions
// for the cube.
size_t TiledSupercell::FindNearestNeighborConfigIndex(
    const size_t referenceLatticeId,    // site in the current small config
    const size_t neighborLatticeId,     // site whose owning config index must be determined
    const size_t configIdxInCube) const // index of the current small config in which reference Site lies
{
  const double cubeSize = static_cast<double>(cubeObj_.GetSizeOfCube());

  // Reference scaled position (inside the cube tiling)
  const Vector3d refOffset = cubeObj_.GetRelativePosition(configIdxInCube).cast<double>();
  const Vector3d refPosition = (smallConfig_.GetRelativePositionOfLattice(referenceLatticeId) + refOffset) / cubeSize;

  double minDistance = std::numeric_limits<double>::max();
  // This is the index of neighbouring small config in which
  // neighborLatticeId lies closest to reference lattice Id
  size_t bestConfigIdx = std::numeric_limits<size_t>::max();

  // Loop over neighbor configs in the cube
  for (size_t nnIdx : cubeObj_.GetNeighbors(configIdxInCube))
  {
    const Vector3d nnOffset = cubeObj_.GetRelativePosition(nnIdx).cast<double>();
    const Vector3d nnPosition = (smallConfig_.GetRelativePositionOfLattice(neighborLatticeId) + nnOffset) / cubeSize;

    // Minimum image convention (PBC)
    Vector3d diff = nnPosition - refPosition;
    for (int i = 0; i < 3; ++i)
      diff[i] -= std::round(diff[i]); // Wrap into [-0.5, 0.5]

    const double distance = diff.norm();
    if (distance < minDistance)
    {
      minDistance = distance;
      bestConfigIdx = nnIdx;
    }
  }

  return bestConfigIdx;
}

vector<pair<size_t, size_t>> TiledSupercell::GetLatticeNeighbors(
    const size_t latticeId,
    const size_t configIdxInCube,
    const size_t maxBondOrder) const
{
  // Neighbour list of a lattice Id in the small config
  auto neighbourLatticeIdsVector = smallConfig_.GetNeighborLatticeIdVectorOfLattice(
      latticeId,
      maxBondOrder);

  // Store <neighbourLatticeId, smallConfigIdx>
  // smallConfigIdx is the index of smallConfig in the cube
  vector<pair<size_t, size_t>> latticeIdNeighbourVector;
  latticeIdNeighbourVector.reserve(neighbourLatticeIdsVector.size());

  // offset of the current smallConfig in the cube
  Vector3d offset = cubeObj_.GetRelativePosition(configIdxInCube).cast<double>();

  // Iterate over the nnLatticeIds

  for (const size_t &nnLatticeId : neighbourLatticeIdsVector)
  {
    // Now check whether this nnLatticeId lies in same config
    bool withInSameConfig = LatticeInSameConfig(
        latticeId,
        nnLatticeId);

    if (!withInSameConfig)
    {
      // If the nnLatticeId does not lies in the same smallConfig_ then get the
      // index of the smallConfig inside the cube in which the nnLatticeId lies
      size_t nnLatticeIdConfigIdx = FindNearestNeighborConfigIndex(
          latticeId,
          nnLatticeId,
          configIdxInCube);

      // <LATTICE_ID_IN_SMALL_CFG , SMALL_CFG_ID_IN_CUBE>
      latticeIdNeighbourVector.emplace_back(
          make_pair(nnLatticeId, nnLatticeIdConfigIdx));
    }
    else
    {
      // If it lies in the same config then assign the same smallCfg Index
      // <LATTICE_ID_IN_SMALL_CFG , SMALL_CFG_ID_IN_CUBE>
      latticeIdNeighbourVector.emplace_back(
          make_pair(nnLatticeId, configIdxInCube));
    }
  }

  // <LATTICE_ID_IN_SMALL_CFG , SMALL_CFG_ID_IN_CUBE>
  return latticeIdNeighbourVector;
}

void TiledSupercell::InitializeAtomVector()
{
  atomVector_.reserve(totalNumOfSites_);

  for (size_t i = 0; i < totalNumOfSites_; i++)
  {
    auto latticeAndConfigIdPair = GetLatticeAndConfigIdFromAtomId(i);

    auto latticeId = latticeAndConfigIdPair.first;
    // auto smallConfigId = configAndLatticeIdPair.second;

    atomVector_.emplace_back(
        smallConfig_.GetElementOfLattice(latticeId));
  }
}

void TiledSupercell::UpdateNeighbourLists(
    const size_t maxBondOrder)
{
  // The reference smallConfig sits at configIndex 0
  // It does matter what index one assign to the reference config
  // because the encoded Index will be used which applies to all the smallConfig
  // inside the cube and one only need to know about the index of current smallConfig
  const size_t referenceConfigIdx = 0;
  // Get the neighbour list of the current small config index in cube
  auto neighbourOfSmallConfig = cubeObj_.GetNeighbors(referenceConfigIdx);

  // Map of the smallConfigIdx to the position index in the sorted neighbour vector
  // for a give smallConfigIdx in the cube
  unordered_map<size_t, int> smallConfigIdxToEncodedIdxMap;
  smallConfigIdxToEncodedIdxMap.reserve(neighbourOfSmallConfig.size());

  for (int idx = 0; idx < neighbourOfSmallConfig.size(); idx++)
  {
    smallConfigIdxToEncodedIdxMap[neighbourOfSmallConfig[idx]] = idx;
  }

  neighbourList_.reserve(maxBondOrder);

  vector<vector<pair<size_t, int>>> neighbourListBondOrder;
  neighbourListBondOrder.reserve(numSitesPerSmallConfig_);

  vector<pair<size_t, int>> encodedNeighbourPairVector;

  for (size_t bondOrder = 1; bondOrder < maxBondOrder + 1; bondOrder++)
  {
    neighbourListBondOrder.clear();

    for (size_t latticeId = 0; latticeId < numSitesPerSmallConfig_; latticeId++)
    {
      // LATTICE_ID_IN_SMALL_CFG : {<NN_LATTICE_ID_IN_SMALL_CFG , SMALL_CFG_ID_IN_CUBE>}
      auto neighbourPairVector = GetLatticeNeighbors(
          latticeId,
          referenceConfigIdx,
          bondOrder);

      // encode the neighbour vector
      encodedNeighbourPairVector.clear();
      encodedNeighbourPairVector.reserve(neighbourPairVector.size());

      for (const auto &entry : neighbourPairVector)
      {
        size_t nnlatticeId = entry.first;
        size_t smallConfigIdxInCube = entry.second;

        if (smallConfigIdxInCube != referenceConfigIdx)
        {
          // Get the encoded index for the map
          int encodedIndex = smallConfigIdxToEncodedIdxMap.at(smallConfigIdxInCube);

          encodedNeighbourPairVector.emplace_back(
              make_pair(nnlatticeId,
                        encodedIndex));
        }
        else
        {
          // If the nnLatticeId for a give latticeId lies in same config then
          encodedNeighbourPairVector.emplace_back(
              make_pair(nnlatticeId,
                        -1));
        }
      }
      neighbourListBondOrder.emplace_back(encodedNeighbourPairVector);
    }
    neighbourList_.emplace_back(neighbourListBondOrder);
  }
}

void TiledSupercell::PrintTiledSupercell() const
{
  const int width = 80;
  const int labelWidth = 30;
  const int valueWidth = width - labelWidth - 2; // 2 for spacing

  // Header
  cout << string(width, '-') << "\n";
  cout << setw((width + 25) / 2) << right << "TiledSupercell Info" << "\n";
  cout << string(width, '-') << "\n";

  // Basic info
  cout << left << setw(labelWidth) << "Num small configs:" << numSmallConfig_ << "\n";
  cout << left << setw(labelWidth) << "Num sites per small config:" << numSitesPerSmallConfig_ << "\n";
  cout << left << setw(labelWidth) << "Total number of sites:" << totalNumOfSites_ << "\n";
  cout << left << setw(labelWidth) << "Atom vector size:" << atomVector_.size() << "\n";

  cout << string(width, '-') << "\n\n";
}
