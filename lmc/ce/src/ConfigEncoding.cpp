#include "ConfigEncoding.h"

static std::unordered_map<ClusterType, size_t, boost::hash<ClusterType>> ConvertSetToHashMap(
    const std::set<ClusterType> &cluster_type_set)
{
  std::unordered_map<ClusterType, size_t, boost::hash<ClusterType>> cluster_type_count;
  for (const auto &cluster_type : cluster_type_set)
  {
    cluster_type_count[cluster_type] = 0;
  }
  return cluster_type_count;
}

ConfigEncoding::ConfigEncoding(
    const Config &config,
    const set<Element> &elementSet,
    const size_t maxBondOrder,
    const size_t maxClusterSize) : maxBondOrder_(maxBondOrder),
                                   maxClusterSize_(maxClusterSize),
                                   initializedClusterHashSet_(
                                       InitializeClusterTypeSet(
                                           config,
                                           elementSet,
                                           maxClusterSize,
                                           maxBondOrder)),
                                   clusterTypeCountHashMap_(
                                       ConvertSetToHashMap(initializedClusterHashSet_)),
                                   latticeClusterCountHashMap_(
                                       CountLatticeClusterTypes(
                                           config,
                                           maxClusterSize,
                                           maxBondOrder)),
                                   allLatticeClusterHashSet_(
                                       FindAllLatticeClusters(
                                           config,
                                           maxClusterSize,
                                           maxBondOrder,
                                           {}))

{
  std::cout << "Cluster Types and Corresponding LatticeClusterCount:\n";
  for (const auto &clusterType : initializedClusterHashSet_)
  {
    std::cout << clusterType << " : "
              << latticeClusterCountHashMap_.at(clusterType.lattice_cluster_type_)
              << "\n";
  }
}

VectorXd ConfigEncoding::GetEncodeVector(const Config &config)
{
  auto clusterTypeCountHashMap = clusterTypeCountHashMap_;
  for (const auto &latticeCluster : allLatticeClusterHashSet_)
  {
    auto atomClusterType = IdentifyAtomClusterType(config, latticeCluster.GetLatticeIdVector());
    clusterTypeCountHashMap.at(ClusterType(atomClusterType, latticeCluster.GetClusterType()))++;
  }

  VectorXd encodeVector(initializedClusterHashSet_.size());
  int idx = 0;
  for (const auto &clusterType : initializedClusterHashSet_)
  {
    auto clusterTypeCount = clusterTypeCountHashMap.at(clusterType);

    // Count of Cluster Types for given configuration
    auto countCluster = static_cast<double>(clusterTypeCount);
    // Count of Cluster Types for normalization
    auto totalLatticeCluster = static_cast<double>(latticeClusterCountHashMap_.at(clusterType.lattice_cluster_type_));

    encodeVector(idx) = countCluster / totalLatticeCluster;

    // std::cout << clusterType << " : " << countCluster << " : " << totalLatticeCluster << std::endl;
    ++idx;
  }

  return encodeVector;
}
