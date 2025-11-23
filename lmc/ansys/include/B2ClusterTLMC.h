#ifndef LMC_ANSYS_INCLUDE_B2CLUSTERTLMC_H_
#define LMC_ANSYS_INCLUDE_B2CLUSTERTLMC_H_

#include <queue>
#include <cstddef>
#include <unordered_set>

#include "Config.h"
#include "UnionFind.h"
#include "TiledSupercell.h"
#include "LatticeSiteMapping.hpp"


using namespace std;

class B2ClusterTLMC
{
public:
  B2ClusterTLMC(
      const TiledSupercell &tiledSupercell);

  void WriteB2ClusterConfig(const string &filename);

  // Returns atomId of the clusters formed
  vector<unordered_set<size_t>> GetB2Clusters();

private:
  void BuildB2Clusters();

  /**
   * @brief Compute all B2 clusters for a given small configuration.
   *
   * This function identifies every B2-centered cluster associated with a
   * specified small configuration. For each cluster, it collects:
   *
   * - **fullB2ClusterVector**:
   *   A list of clusters where each entry contains all lattice sites
   *   belonging to the cluster, including neighbors up to the 2nd nearest
   *   neighbor.
   *
   * - **b2CenterClusterVector**:
   *   A list of clusters where each entry contains only the B2 center sites
   *   corresponding to the clusters in @p fullB2ClusterVector.
   *   These reduced clusters are used to merge clusters later and avoid
   *   double counting.
   *
   * @param smallConfigId
   *        ID of the small configuration for which B2-centered clusters
   *        are to be generated.
   *
   * @param fullB2ClusterVector
   *        Output vector where each element is the complete set of lattice
   *        sites (including 1st and 2nd nearest neighbors) forming a B2 cluster.
   *
   * @param b2CenterClusterVector
   *        Output vector containing only the B2 center sites corresponding to
   *        each full cluster, used for cluster merging and uniqueness checks.
   */
  void GetB2ClustersForSmallConfig(
      const size_t smallConfigId,
      vector<unordered_set<size_t>> &fullB2ClusterVector,
      vector<unordered_set<size_t>> &b2CenterClusterVector);

  bool isB2(const LatticeSiteMapping &latticeSiteId) const;

  const TiledSupercell &tiledSupercell_;
  vector<unordered_set<size_t>> b2ClusterVector_{};
};

#endif // LMC_ANSYS_INCLUDE_B2CLUSTERTLMC_H_