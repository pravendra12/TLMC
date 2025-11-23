#include "B2ClusterTLMC.h"

B2ClusterTLMC::B2ClusterTLMC(const TiledSupercell &tiledSupercell) : tiledSupercell_(tiledSupercell)
{
  BuildB2Clusters();
}

void B2ClusterTLMC::WriteB2ClusterConfig(const string &filename)
{
  auto largeConfig = tiledSupercell_.MakeSupercell();
  auto numSites = tiledSupercell_.GetTotalNumOfSites();

  // Auxiliary lists
  Config::VectorVariant clusterIdVec = vector<int>(numSites, -1);
  Config::VectorVariant clusterSizeVec = vector<size_t>(numSites, 0);
  Config::VectorVariant isSharedVec = vector<int>(numSites, 0); // 0 = not shared, 1 = shared

  unordered_map<size_t, size_t> atomClusterCount;

  // First pass: count how many clusters each atom belongs to
  for (size_t clusterId = 0; clusterId < b2ClusterVector_.size(); ++clusterId)
  {
    for (size_t atomId : b2ClusterVector_[clusterId])
    {
      atomClusterCount[atomId]++;
    }
  }

  // Second pass: fill auxiliary lists
  for (size_t clusterId = 0; clusterId < b2ClusterVector_.size(); ++clusterId)
  {
    size_t clusterSize = b2ClusterVector_[clusterId].size();
    for (size_t atomId : b2ClusterVector_[clusterId])
    {
      get<vector<int>>(clusterIdVec)[atomId] = int(clusterId);
      get<vector<size_t>>(clusterSizeVec)[atomId] = clusterSize;

      if (atomClusterCount[atomId] > 1)
        get<vector<int>>(isSharedVec)[atomId] = 1; // mark as shared
    }
  }

  map<string, Config::VectorVariant> auxiliaryLists;
  auxiliaryLists["clusterId"] = clusterIdVec;
  auxiliaryLists["clusterSize"] = clusterSizeVec;
  auxiliaryLists["sharedAtom"] = isSharedVec; // new shared flag

  map<string, Config::ValueVariant> globalList; // empty

  Config::WriteXyzExtended(filename, largeConfig, auxiliaryLists, globalList);
}

vector<unordered_set<size_t>> B2ClusterTLMC::GetB2Clusters()
{
  return b2ClusterVector_;
}

void B2ClusterTLMC::BuildB2Clusters()
{
  // Full clusters (centers + 1NN + 2NN)
  vector<vector<unordered_set<size_t>>> fullClustersPerTile;

  // Contains ONLY the true B2 centers, used for merging logic
  vector<vector<unordered_set<size_t>>> b2CentersPerTile;

  const auto numOfSmallConfigs = tiledSupercell_.GetNumOfSmallConfig();

  // Build the clusters tile-by-tile
  /*
  // Series version
  for (size_t smallConfigId = 0; smallConfigId < numOfSmallConfigs; smallConfigId++)
  {
    vector<unordered_set<size_t>> fullB2ClusterVector;
    vector<unordered_set<size_t>> b2CenterClusterVector;

    // Modified GetB2ClustersForSmallConfig fills BOTH outputs
    GetB2ClustersForSmallConfig(
        smallConfigId,
        fullB2ClusterVector,
        b2CenterClusterVector);

    fullClustersPerTile.emplace_back(std::move(fullB2ClusterVector));
    b2CentersPerTile.emplace_back(std::move(b2CenterClusterVector));
  }
  */

  fullClustersPerTile.resize(numOfSmallConfigs);
  b2CentersPerTile.resize(numOfSmallConfigs);

#pragma omp parallel for
  for (size_t smallConfigId = 0; smallConfigId < numOfSmallConfigs; smallConfigId++)
  {
    vector<unordered_set<size_t>> fullB2ClusterVector;
    vector<unordered_set<size_t>> b2CenterClusterVector;

    GetB2ClustersForSmallConfig(
        smallConfigId,
        fullB2ClusterVector,
        b2CenterClusterVector);

    fullClustersPerTile[smallConfigId] = std::move(fullB2ClusterVector);
    b2CentersPerTile[smallConfigId] = std::move(b2CenterClusterVector);
  }

  // Union–Find ONLY on B2 CENTERS (NOT neighbors!)

  UnionFind<size_t> uf;

  // Each cluster index corresponds to the same B2-center group
  for (size_t tile = 0; tile < numOfSmallConfigs; tile++)
  {
    const auto &b2CenterClusters = b2CentersPerTile[tile];

    for (const auto &b2CenterCluster : b2CenterClusters)
    {
      if (b2CenterCluster.empty())
        continue;

      // Pick the first B2 center as the root
      auto it = b2CenterCluster.begin();
      size_t root = *it++;

      // Union only the B2 centers
      for (; it != b2CenterCluster.end(); ++it)
        uf.union_sets(root, *it);
    }
  }

  // Build merged clusters: union-find defines connectivity,
  // but we include ALL atoms from the full clusters

  unordered_map<size_t, unordered_set<size_t>> mergedClusters;

  for (size_t tile = 0; tile < numOfSmallConfigs; tile++)
  {
    const auto &fullClusters = fullClustersPerTile[tile];
    const auto &b2CenterClusters = b2CentersPerTile[tile];

    for (size_t i = 0; i < fullClusters.size(); i++)
    {
      const auto &fullCluster = fullClusters[i];
      const auto &centerCluster = b2CenterClusters[i];

      if (centerCluster.empty())
        continue;

      // All B2 centers in this cluster share the same UF root
      size_t rootCenter = uf.find(*centerCluster.begin());

      // Insert ALL atoms of the full cluster into final merged cluster
      for (size_t atomId : fullCluster)
        mergedClusters[rootCenter].insert(atomId);
    }
  }

  // Write into the class storage (final B2 clusters)
  b2ClusterVector_.clear();
  for (auto &entry : mergedClusters)
    b2ClusterVector_.emplace_back(std::move(entry.second));
}

// This function will return the B2 cluster for a given small config Id
void B2ClusterTLMC::GetB2ClustersForSmallConfig(
    const size_t smallConfigId,
    vector<unordered_set<size_t>> &fullB2ClusterVector,
    vector<unordered_set<size_t>> &b2CenterClusterVector)
{

  // These will store all the sites in cluster including its neighbours (upto 2nd nn)
  // vector<unordered_set<size_t>> fullB2ClusterVector;

  // Will contain only the B2 centers corresponding to the above clusters which
  // will be used to merge the clusters later to avoid double counting
  // vector<unordered_set<size_t>> b2CenterClusterVector;

  const size_t numSitesInSmallConfig = tiledSupercell_.GetNumOfSitesPerSmallConfig();

  unordered_set<size_t> visitedB2; // only B2 sites marked visited

  // latticeId still starts from 0 but the LatticeSiteMapping will be different
  // for different smallConfigId
  for (size_t latticeId = 0; latticeId < numSitesInSmallConfig; ++latticeId)
  {
    const LatticeSiteMapping latticeSiteId(latticeId, smallConfigId);

    // atom Id is unique so it will be used for keeping the track whether a site
    // has been visited or not Also the atomId will be used for clusters as well
    const auto atomId = tiledSupercell_.GetAtomIdFromLatticeAndConfigId(latticeSiteId);

    // skip if site already processed as B2 center
    if (visitedB2.count(atomId))
      continue;

    // check if this site is a B2 center
    if (!isB2(latticeSiteId))
      continue;

    // start new cluster
    unordered_set<size_t> currentCluster;

    // Will only store the b2 centers
    unordered_set<size_t> b2Centers;

    queue<size_t> bfsQueue;

    bfsQueue.push(atomId);
    visitedB2.insert(atomId);

    while (!bfsQueue.empty())
    {
      size_t center = bfsQueue.front();
      bfsQueue.pop();

      // add center to cluster
      currentCluster.insert(center);
      b2Centers.insert(center);

      // get 1NN + 2NN neighbors
      // now go back to LatticeSiteMapping from atomId and then get the neighbours

      auto centerLatticeSiteId = tiledSupercell_.GetLatticeSiteMappingFromAtomId(center);

      // get neighbours of the center's small config
      const auto &neighboursOfCenter = tiledSupercell_.GetCube().GetNeighbors(centerLatticeSiteId.smallConfigId);

      auto nnEncodedUpToSecond = tiledSupercell_.GetNeighborLatticeIdsUpToOrder(
          centerLatticeSiteId.latticeId, 2);

      // add all NN (B2 or not) into the cluster
      for (const auto &nnEncodedId : nnEncodedUpToSecond)
      {
        size_t nnSmallConfigId = centerLatticeSiteId.smallConfigId;

        if (nnEncodedId.encodedSmallConfigId != -1)
        {
          nnSmallConfigId = neighboursOfCenter[size_t(nnEncodedId.encodedSmallConfigId)];
        }

        const LatticeSiteMapping nnLatticeSiteId(nnEncodedId.latticeId, nnSmallConfigId);
        const size_t nnAtomId = tiledSupercell_.GetAtomIdFromLatticeAndConfigId(nnLatticeSiteId);

        currentCluster.insert(nnAtomId);

        if (visitedB2.count(nnAtomId))
          continue;

        if (isB2(nnLatticeSiteId))
        {
          visitedB2.insert(nnAtomId);
          bfsQueue.push(nnAtomId);
        }
      }
    }

    fullB2ClusterVector.push_back(std::move(currentCluster));
    b2CenterClusterVector.push_back(std::move(b2Centers));
  }
}

bool B2ClusterTLMC::isB2(const LatticeSiteMapping &latticeSiteId) const
{
  const Element centralElement =
      tiledSupercell_.GetElementAtSite(latticeSiteId);

  if (centralElement == Element("X"))
    return false;

  // Neighbours of the current cube
  const vector<size_t> &neighboursOfSmallConfig =
      tiledSupercell_.GetCube().GetNeighbors(latticeSiteId.smallConfigId);

  //  These are encoding
  const auto firstNN = tiledSupercell_.GetNeighborLatticeIdVectorOfLattice(
      latticeSiteId.latticeId, 1);

  const auto latticeSiteEncodingNN1 = firstNN[0];
  size_t nnSmallConfigId = latticeSiteId.smallConfigId;

  if (latticeSiteEncodingNN1.encodedSmallConfigId != -1)
  {
    nnSmallConfigId = neighboursOfSmallConfig[size_t(latticeSiteEncodingNN1.encodedSmallConfigId)];
  }

  // The 1NN element type must be uniform and opposite to central
  const Element firstNNElement = tiledSupercell_.GetElementAtSite(
      LatticeSiteMapping(latticeSiteEncodingNN1.latticeId,
                         nnSmallConfigId));

  if (firstNNElement == centralElement)
    return false;

  // First NN conditions
  for (const auto &nnLatticeSiteEncoding : firstNN)
  {
    nnSmallConfigId = latticeSiteId.smallConfigId;

    if (nnLatticeSiteEncoding.encodedSmallConfigId != -1)
    {
      nnSmallConfigId = neighboursOfSmallConfig[size_t(nnLatticeSiteEncoding.encodedSmallConfigId)];
    }

    auto nnElement = tiledSupercell_.GetElementAtSite(LatticeSiteMapping(
        nnLatticeSiteEncoding.latticeId,
        nnSmallConfigId));

    if (nnElement != firstNNElement)
      return false;
  }

  const auto secondNN = tiledSupercell_.GetNeighborLatticeIdVectorOfLattice(
      latticeSiteId.latticeId, 2);

  // Second NN conditions
  for (const auto &nnLatticeSiteEncoding : secondNN)
  {
    nnSmallConfigId = latticeSiteId.smallConfigId;

    if (nnLatticeSiteEncoding.encodedSmallConfigId != -1)
    {
      nnSmallConfigId =
          neighboursOfSmallConfig[size_t(nnLatticeSiteEncoding.encodedSmallConfigId)];
    }

    auto nnElement = tiledSupercell_.GetElementAtSite(LatticeSiteMapping(
        nnLatticeSiteEncoding.latticeId,
        nnSmallConfigId));

    if (nnElement != centralElement)
      return false;
  }

  return true;
}
