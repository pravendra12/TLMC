#ifndef LMC_CE_INCLUDE_CONFIGENCODING_H_
#define LMC_CE_INCLUDE_CONFIGENCODING_H_

#include <Eigen/Dense>
#include "Config.h"
#include "ClusterExpansion.h"

using namespace std;
using namespace Eigen;

class ConfigEncoding
{
  public:
  ConfigEncoding(
      const Config &config,
      const set<Element> &elementSet,
      const size_t maxBondOrder,
      const size_t maxClusterSize);

  VectorXd GetEncodeVector(const Config &config);


  private:
    const size_t maxBondOrder_;
    const size_t maxClusterSize_;

    const set<Element> elementSet_;

    const set<ClusterType> initializedClusterHashSet_;

    const unordered_map<ClusterType, size_t, boost::hash<ClusterType>> clusterTypeCountHashMap_;

    const unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>> latticeClusterCountHashMap_;

    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> allLatticeClusterHashSet_;

  
};

#endif // LMC_CE_INCLUDE_CONFIGENCODING_H_
