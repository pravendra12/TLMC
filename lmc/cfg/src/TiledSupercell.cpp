#include "TiledSupercell.h"
#include "Config.h"

TiledSupercell::TiledSupercell(
    Config &smallConfig,
    const Cube &cubeObj) : smallConfig_(smallConfig),
                           cubeObj_(move(cubeObj)),
                           superBasis_(smallConfig.GetBasis() * static_cast<double>(cubeObj_.GetSizeOfCube())),
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

const Config &TiledSupercell::GetSmallConfig() const
{
  return smallConfig_;
}

const Matrix3d &TiledSupercell::GetSuperBasis() const
{
  return superBasis_;
}

const Cube &TiledSupercell::GetCube() const
{
  return cubeObj_;
}

size_t TiledSupercell::GetNumOfSitesPerSmallConfig() const
{
  return numSitesPerSmallConfig_;
}

size_t TiledSupercell::GetNumOfSmallConfig() const
{
  return numSmallConfig_;
}

size_t TiledSupercell::GetTotalNumOfSites() const
{
  return totalNumOfSites_;
}

// these latticePair are bascially teh
// <LatticeId, smallConfigId> which will be used to
// get the atomId then the swap will happen
void TiledSupercell::LatticeJump(
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteMappingPair)
{
  auto latticeSiteMapping1 = latticeSiteMappingPair.first;
  auto latticeSiteMapping2 = latticeSiteMappingPair.second;

  auto atomId1 = GetAtomIdFromLatticeAndConfigId(
      latticeSiteMapping1);
  auto atomId2 = GetAtomIdFromLatticeAndConfigId(
      latticeSiteMapping2);

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
  // Also make a swap to atomIndicesVector_
  std::swap(atomIndexVector_[atomId1], atomIndexVector_[atomId2]);
}

vector<LatticeSiteEncodedMapping> TiledSupercell::GetNeighborLatticeIdVectorOfLattice(
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

vector<LatticeSiteEncodedMapping> TiledSupercell::GetNeighborLatticeIdsUpToOrder(
    const size_t &latticeId,
    const size_t &maxBondOrder) const
{
  if (maxBondOrder == 0)
  {
    throw std::invalid_argument(
        "Error in TiledSupercell::GetNeighborLatticeIdsUpToOrder: "
        "maxBondOrder must be >= 1");
  }
  if (maxBondOrder > neighbourList_.size())
  {
    throw std::out_of_range(
        "Error in TiledSupercell::GetNeighborLatticeIdsUpToOrder: "
        "maxBondOrder (" +
        std::to_string(maxBondOrder) +
        ") exceeds neighbourList_ size (" +
        std::to_string(neighbourList_.size()) + ")");
  }
  if (latticeId >= neighbourList_[0].size())
  {
    throw std::out_of_range(
        "Error in TiledSupercell::GetNeighborLatticeIdsUpToOrder: "
        "latticeId (" +
        std::to_string(latticeId) +
        ") exceeds number of lattice sites (" +
        std::to_string(neighbourList_[0].size()) + ")");
  }

  vector<LatticeSiteEncodedMapping> neighborsUpToOrder;
  for (size_t order = 1; order < maxBondOrder + 1; ++order)
  {
    const auto &ids = GetNeighborLatticeIdVectorOfLattice(latticeId, order);
    neighborsUpToOrder.insert(neighborsUpToOrder.end(), ids.begin(), ids.end());
  }

  return neighborsUpToOrder;
}

vector<NeighbourOfPair> TiledSupercell::GetNeighboringLatticeIdSetOfPair(
    const pair<size_t, size_t> &latticeIdPair,
    const size_t &maxBondOrder,
    const bool removeLatticeIdDuplicates) const
{

  // MAJOR ISSUE

  // 1
  // common neighbours are still a issue need to think about the common neibhours
  // as two same latticeId can have diffent encodedIndex it like seeing the same id
  // from two different sides, this will be a major issue, but in principle these
  // common sites viewed from either site will have same element so keep one of it

  // SOLUTION
  // This is solved using NeighbourOfPair which also store this neighbour belongs to
  // which site in the latticeIdPair first or second

  // 2
  // Right now this function is only be expected to be used if the first nn latticeIdPair in the
  // TileSupercell, for the case where they are not the first nn in the TiledSupercell
  // then return all the neighbours without any use seen set
  // One way to make it more generic to use a bool removeLatticeIdDuplicates

  // SOLUTION
  // Use removeLatticeIdDuplicates tag which will specify whether to remove those
  // sites which have same latticeId but different encodedConfigIdx, it is useful
  // when one have to get the neighbours of pair which contains latticeIds at first nn
  // in TiledSupercell, also there can be case when the two same latticeId may not be
  // first nn in the TiledSupercell in that case one would be interested to get all the
  // neighbours, not sure whether later will be used at all but a now the function is
  // more generic

  vector<NeighbourOfPair> neighbours;

  // If deduplication is needed, use a set to track added latticeIds
  unordered_set<size_t> seen;

  for (size_t bondOrder = 1; bondOrder <= maxBondOrder; ++bondOrder)
  {
    auto neighborsVector1 = GetNeighborLatticeIdVectorOfLattice(latticeIdPair.first, bondOrder);
    auto neighborsVector2 = GetNeighborLatticeIdVectorOfLattice(latticeIdPair.second, bondOrder);

    for (const auto &nn : neighborsVector1)
    {
      auto id = nn.latticeId;
      if (id != latticeIdPair.first && id != latticeIdPair.second)
      {
        if (!removeLatticeIdDuplicates || seen.insert(id).second)
          neighbours.emplace_back(nn, NeighbourOfPair::SourceSite::First);
      }
    }

    for (const auto &nn : neighborsVector2)
    {
      auto id = nn.latticeId;
      if (id != latticeIdPair.first && id != latticeIdPair.second)
      {
        if (!removeLatticeIdDuplicates || seen.insert(id).second)
          neighbours.emplace_back(nn, NeighbourOfPair::SourceSite::Second);
      }
    }
  }

  return neighbours;
}

void TiledSupercell::SetElementAtSite(
    const LatticeSiteMapping &latticeSiteMapping,
    const Element &element)
{
  size_t atomId = GetAtomIdFromLatticeAndConfigId(
      latticeSiteMapping);

  if (atomId >= totalNumOfSites_)
  {
    throw std::out_of_range(
        "TiledSupercell::SetElementAtSite error: "
        "atomId " +
        std::to_string(atomId) +
        " (from latticeId " + std::to_string(latticeSiteMapping.latticeId) +
        ", smallConfigId " + std::to_string(latticeSiteMapping.smallConfigId) +
        ") exceeds total number of sites " + std::to_string(totalNumOfSites_));
  }
  atomVector_[atomId] = element;
  atomIndexVector_[atomId] = static_cast<uint64_t>(element.GetAtomicIndex());
}

Element TiledSupercell::GetElementAtSite(
    const LatticeSiteMapping &latticeSiteMapping) const
{
  size_t atomId = GetAtomIdFromLatticeAndConfigId(
      latticeSiteMapping);

  if (atomId >= totalNumOfSites_)
  {
    throw std::out_of_range(
        "TiledSupercell::GetElementAtSite error: "
        "atomId " +
        std::to_string(atomId) +
        " (from latticeId " + std::to_string(latticeSiteMapping.latticeId) +
        ", smallConfigId " + std::to_string(latticeSiteMapping.smallConfigId) +
        ") exceeds total number of sites " + std::to_string(totalNumOfSites_));
  }

  return atomVector_.at(atomId);
}

LatticeSiteMapping TiledSupercell::GetVacancySiteId() const
{
  size_t xCount = std::count(atomVector_.begin(), atomVector_.end(), Element("X"));

  if (xCount == 0)
  {
    throw std::runtime_error(
        "Error in `TiledSupercell::GetVacancySiteId`: No vacancy found in TiledSupercell. Expected exactly one vacancy.");
  }

  if (xCount > 1)
  {
    throw std::runtime_error(
        "Error in `TiledSupercell::GetVacancySiteId`: More than one vacancy found in TiledSupercell. Expected exactly one vacancy.");
  }

  // Safe to get the index now
  auto it = std::find(atomVector_.begin(), atomVector_.end(), Element("X"));
  size_t vacancyIndex = std::distance(atomVector_.begin(), it);

  auto vacancySiteId = GetLatticeSiteMappingFromAtomId(vacancyIndex);

  return vacancySiteId;
}

LatticeSiteMapping TiledSupercell::GetLatticeSiteMappingFromEncoding(
    const LatticeSiteEncodedMapping &latticeSiteEncoding,
    const LatticeSiteMapping &refLatticeSiteMapping) const
{
  size_t neighbourConfigIdx;
  if (latticeSiteEncoding.encodedSmallConfigId == -1)
  {
    neighbourConfigIdx = refLatticeSiteMapping.smallConfigId;
  }
  else
  {
    const auto &neighborsOfSmallConfig = cubeObj_.GetNeighbors(refLatticeSiteMapping.smallConfigId);

    neighbourConfigIdx = neighborsOfSmallConfig[latticeSiteEncoding.encodedSmallConfigId];
  }

  return LatticeSiteMapping(
      latticeSiteEncoding.latticeId,
      neighbourConfigIdx);
}

// Return the absolute relative position vector of a lattice site in the supercell
Vector3d TiledSupercell::GetRelativePositionOfLatticeSiteMapping(
    const LatticeSiteMapping &latticeSiteId) const
{
  const size_t smallConfigIdx = latticeSiteId.smallConfigId;
  const size_t latticeId = latticeSiteId.latticeId;

  // Offset of this smallConfig in the cube
  Vector3d offset = cubeObj_.GetRelativePosition(smallConfigIdx).cast<double>();

  // Position of latticeId in the small config
  Vector3d smallConfigPos = smallConfig_.GetRelativePositionMatrix().col(static_cast<int>(latticeId));

  // Supercell relative position = small config position + offset (scaled by cube size)
  Vector3d relPos = (smallConfigPos + offset) / static_cast<double>(cubeObj_.GetSizeOfCube());

  return relPos;
}

// Return the relative distance vector between two lattice sites in the supercell
Vector3d TiledSupercell::GetRelativeDistanceVectorLattice(
    const LatticeSiteMapping &site1,
    const LatticeSiteMapping &site2) const
{
  Vector3d relVec = GetRelativePositionOfLatticeSiteMapping(site2) -
                    GetRelativePositionOfLatticeSiteMapping(site1);

  // Apply periodic boundary conditions
  for (int kDim = 0; kDim < 3; ++kDim)
  {
    while (relVec[kDim] >= 0.5)
      relVec[kDim] -= 1.0;
    while (relVec[kDim] < -0.5)
      relVec[kDim] += 1.0;
  }
  return relVec;
}

LatticeSiteMapping TiledSupercell::GetLatticeSiteMappingFromAtomId(
    const size_t &atomId) const
{
  if (atomId >= totalNumOfSites_)
  {
    throw std::out_of_range("Error in `TiledSupercell::GetLatticeSiteMappingFromAtomId`: Atom index exceeds total number of atoms");
  }

  size_t smallConfigId = atomId / numSitesPerSmallConfig_;
  size_t latticeId = atomId % numSitesPerSmallConfig_;

  return LatticeSiteMapping{latticeId, smallConfigId};
}

size_t TiledSupercell::GetAtomIdFromLatticeAndConfigId(
    const LatticeSiteMapping &latticeSiteMapping) const
{

  size_t atomId = latticeSiteMapping.smallConfigId * numSitesPerSmallConfig_ + latticeSiteMapping.latticeId;

  if (atomId >= totalNumOfSites_)
  {
    throw std::out_of_range("Error in `TiledSupercell::GetAtomIdFromLatticeAndConfigId` : Computed atom index exceeds total number of atoms");
  }
  return atomId;
}

struct Vector3iHash
{
  size_t operator()(const Eigen::Vector3i &v) const
  {
    return std::hash<int>()(v[0]) ^ std::hash<int>()(v[1] << 1) ^ std::hash<int>()(v[2] << 2);
  }
};

void TiledSupercell::UpdateAtomVector(
    const Config &config)
{
  if (config.GetNumAtoms() != totalNumOfSites_)
    throw std::runtime_error("Error in `TiledSupercell::UpdateAtomVector`: Config size does not match Tiled Supercell size");

  const auto &relPos = config.GetRelativePositionMatrix(); // 3 x numAtomsSmallCfg

  // Build a hash map: discretized positions → lattice index
  std::unordered_map<Eigen::Vector3i, size_t, Vector3iHash> posToLatticeId;
  for (size_t j = 0; j < config.GetNumLattices(); ++j)
  {
    Eigen::Vector3i key = (relPos.col(j) / constants::kEpsilon).array().round().cast<int>();
    posToLatticeId[key] = j;
  }

  for (size_t i = 0; i < totalNumOfSites_; ++i)
  {
    auto latticeSiteMapping = GetLatticeSiteMappingFromAtomId(i);
    Vector3d pos = GetRelativePositionOfLatticeSiteMapping(latticeSiteMapping);

    Eigen::Vector3i key = (pos / constants::kEpsilon).array().round().cast<int>();

    auto it = posToLatticeId.find(key);
    if (it == posToLatticeId.end())
      throw std::runtime_error("No matching lattice site found");

    SetElementAtSite(latticeSiteMapping, config.GetElementOfLattice(it->second));
  }
}

void TiledSupercell::UpdateAtomVector(
    const vector<uint64_t> &atomicIndicesVector)
{
  if (totalNumOfSites_ != atomicIndicesVector.size())
    throw std::runtime_error("Error in `TiledSupercell::UpdateAtomVector`: Size of atomicIndicesVector does not matches totalNumOfSites in TiledSupercell");

  atomVector_.clear();
  atomIndexVector_.clear();

  atomVector_.reserve(totalNumOfSites_);
  atomIndexVector_.reserve(totalNumOfSites_);

  for (size_t  i = 0; i < totalNumOfSites_; i++)
  {
    auto atomicIndex = size_t(atomicIndicesVector[i]); // AtomicNumber
    atomVector_.emplace_back(Element(atomicIndex));
    atomIndexVector_.emplace_back(atomicIndicesVector[i]);
  }
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
    auto latticeAndConfigIdPair = GetLatticeSiteMappingFromAtomId(i);

    size_t latticeId = latticeAndConfigIdPair.latticeId;
    size_t smallCfgIdx = latticeAndConfigIdPair.smallConfigId;

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
    auto latticeAndConfigIdPair = GetLatticeSiteMappingFromAtomId(i);

    size_t latticeId = latticeAndConfigIdPair.latticeId;
    size_t smallCfgIdx = latticeAndConfigIdPair.smallConfigId;

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

void TiledSupercell::WriteAtomicIndicesToFile(
    const string &filename) const
{

  std::ofstream ofs(filename, std::ios_base::out | std::ios_base::binary);
  if (!ofs)
    throw std::runtime_error("Cannot open file");

  io::filtering_ostream fos;

  // Push compressor if needed
  std::string ext = boost::filesystem::path(filename).extension().string();
  if (ext == ".gz")
  {
    fos.push(io::gzip_compressor());
  }
  else if (ext == ".bz2")
  {
    fos.push(io::bzip2_compressor());
  }

  // Push file stream last
  fos.push(ofs);

  uint64_t numAtoms = atomIndexVector_.size();
  fos.write(reinterpret_cast<const char *>(&numAtoms), sizeof(numAtoms));
  fos.write(reinterpret_cast<const char *>(atomIndexVector_.data()), atomIndexVector_.size() * sizeof(uint64_t));

  // Important: flush to ensure all compressed data is written
  fos.flush();
}

vector<uint64_t> TiledSupercell::ReadAtomicIndicesFromFile(const string &filename)
{

  std::ifstream ifs(filename, std::ios_base::in | std::ios_base::binary);
  if (!ifs)
    throw std::runtime_error("Error in `TiledSupercell::ReadAtomicIndicesFromFile`: Cannot open file for reading");

  io::filtering_istream fis;

  // Push decompressor if needed
  std::string ext = boost::filesystem::path(filename).extension().string();
  if (ext == ".gz")
  {
    fis.push(io::gzip_decompressor());
  }
  else if (ext == ".bz2")
  {
    fis.push(io::bzip2_decompressor());
  }

  // Push file stream last
  fis.push(ifs);

  // Read number of atoms
  uint64_t numAtoms;
  fis.read(reinterpret_cast<char *>(&numAtoms), sizeof(numAtoms));
  if (!fis)
    throw std::runtime_error("Failed to read number of atoms");

  // Allocate vector and read indices
  std::vector<uint64_t> atomicIndexVector(numAtoms);
  fis.read(reinterpret_cast<char *>(atomicIndexVector.data()), numAtoms * sizeof(uint64_t));
  if (!fis)
    throw std::runtime_error("Failed to read atom indices");

  return atomicIndexVector;
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
  atomIndexVector_.reserve(totalNumOfSites_);

  for (size_t i = 0; i < totalNumOfSites_; i++)
  {
    auto latticeAndConfigIdPair = GetLatticeSiteMappingFromAtomId(i);

    auto latticeId = latticeAndConfigIdPair.latticeId;
    // auto smallConfigId = configAndLatticeIdPair.second;

    auto element = smallConfig_.GetElementOfLattice(latticeId);

    atomVector_.emplace_back(
        element);

    atomIndexVector_.emplace_back(
        static_cast<uint64_t>(element.GetAtomicIndex()));
  }
}

void TiledSupercell::UpdateNeighbourLists(
    const size_t maxBondOrder)
{
  neighbourList_.clear();
  neighbourList_.reserve(maxBondOrder);
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

  vector<vector<LatticeSiteEncodedMapping>> neighbourListBondOrder;
  neighbourListBondOrder.reserve(numSitesPerSmallConfig_);

  vector<LatticeSiteEncodedMapping> encodedNeighbourPairVector;

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
        size_t nnLatticeId = entry.first;
        size_t smallConfigIdxInCube = entry.second;

        if (smallConfigIdxInCube != referenceConfigIdx)
        {
          // Get the encoded index for the map
          int encodedIndex = smallConfigIdxToEncodedIdxMap.at(smallConfigIdxInCube);

          encodedNeighbourPairVector.emplace_back(LatticeSiteEncodedMapping(
              nnLatticeId,
              encodedIndex));
        }
        else
        {
          // If the nnLatticeId for a give latticeId lies in same config then
          encodedNeighbourPairVector.emplace_back(
              LatticeSiteEncodedMapping(
                  nnLatticeId,
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
