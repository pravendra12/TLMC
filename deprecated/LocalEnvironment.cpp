#include "LocalEnvironment.h"
#include "Config.h"

Config GetLocalEnvironmentConfig(const Config &config,
                                 size_t latticeId,
                                 const vector<double> &cutoffs)
{
  // config.UpdateNeighborList(cutoffs);

  size_t maxBondOrder = cutoffs.size();
  auto neighbours = config.GetNeighborLatticeIdsUpToOrder(latticeId, maxBondOrder);
  neighbours.emplace_back(latticeId); // include the central atom

  Matrix3d basis = config.GetBasis();
  Matrix3Xd relativePositionMatrix(3, neighbours.size()); // correct resize

  vector<Element> atomVector;
  atomVector.reserve(neighbours.size());

  for (size_t i = 0; i < neighbours.size(); i++)
  {
    size_t id = neighbours[i];
    atomVector.emplace_back(config.GetElementOfLattice(id));

    Vector3d relativePosition = config.GetRelativePositionOfLattice(id);
    relativePositionMatrix.col(i) = relativePosition;
  }

  Config localConfig(basis, relativePositionMatrix, atomVector);

  localConfig.UpdateNeighborList(cutoffs);

  // The config can be centered around (0.5, 0.5, 0.5) if needed
  return localConfig;
}

LocalEnvironment::LocalEnvironment(
    const Config &config,
    size_t latticeId,
    const vector<double> &cutoffs) : config_(GetLocalEnvironmentConfig(config, // copy
                                                                       latticeId,
                                                                       cutoffs))
{}

Config LocalEnvironment::GetLocalConfig()
{

  // size_t vacancyId = config_.GetVacancyLatticeId();
  return config_;

  // Outside the local environment

  //  cout << "---- Inside LE ----" << endl;
  //
  //  cout << "Lattice ID: " << vacancyId << endl;
  //  cout << "Position: " << config_.GetRelativePositionOfLattice(vacancyId).transpose() << endl;
  //  cout << "Element: " << config_.GetElementOfLattice(vacancyId).GetElementString() << endl;
}

VectorXd LocalEnvironment::GetLocalConfigEncoding(
    const size_t maxClusterSize,
    const size_t maxBondOrder)
{

  auto atomVector = config_.GetAtomVector();

  const set<Element> elementSet(atomVector.begin(), atomVector.end());

  auto clusterTypeSet = InitializeClusterTypeSet(
      config_,
      elementSet,
      maxClusterSize,
      maxBondOrder);

  unordered_map<ClusterType, size_t, boost::hash<ClusterType>> clusterTypeCountMap;
  for (const auto &clusterType : clusterTypeSet)
  {
    clusterTypeCountMap[clusterType] = 0;
  }

  auto allLatticeClusterHashSet = FindAllLatticeClusters(
      config_,
      maxClusterSize,
      maxBondOrder, {});

  for (const auto &latticeCluster : allLatticeClusterHashSet)
  {
    auto atomClusterType = IdentifyAtomClusterType(config_, latticeCluster.GetLatticeIdVector());
    clusterTypeCountMap.at(ClusterType(atomClusterType, latticeCluster.GetClusterType()))++;
  }

  VectorXd encodeVector(clusterTypeSet.size());

  auto latticeClusterTypeCount = CountLatticeClusterTypes(
    config_, 
    maxClusterSize, 
    maxBondOrder
  );
  
  int idx = 0;
  for (const auto &clusterType : clusterTypeSet)
  {
    // Count of Cluster Types for given configuration
    auto clusterCount = static_cast<double>(clusterTypeCountMap.at(clusterType));
    auto totalCount = static_cast<double>(latticeClusterTypeCount.at(clusterType.lattice_cluster_type_));

    encodeVector(idx) = clusterCount/totalCount;
    ++idx;

    // cout << clusterType << " : " << clusterCount << " " << totalCount << endl;
  }

  return encodeVector;
}
