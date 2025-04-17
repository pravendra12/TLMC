#include "ClusterDynamics.h"

ClusterDynamics::ClusterDynamics(const Config &config) : config_(config),
                                                         clusters_{},
                                                         visited_{}
{}

void ClusterDynamics::detectB2Clusters(map<string, Config::VectorVariant> &auxiliaryList,
                                       vector<int> &clusterSizeVector)
{
  size_t numAtoms = config_.GetNumAtoms();
  for (size_t i = 0; i < numAtoms; ++i)
  {
    if (visited_.count(i))
      continue;

    unordered_set<size_t> cluster;
    if (growB2Cluster(i, cluster) && !cluster.empty())
    {
      clusters_.push_back(move(cluster));
    }
  }

  mergeAllClusters();

  size_t vacancyId = static_cast<size_t>(-1);

  try
  {
    vacancyId = config_.GetVacancyAtomId();
  }
  catch (const std::exception &e)
  {
    std::cerr << "Warning: No vacancy found, setting vacancyId = -1.\n";
  }

  vector<int> atomClusterTypeVector(numAtoms, -1);
  vector<int> atomClusterSizeVector(numAtoms, -1);

  int clusterId = 1;

  for (const auto &cluster : clusters_)
  {
    int clusterSize = cluster.size();

    // If vacancy is part of the cluster, we exclude it from the count
    if (cluster.find(vacancyId) != cluster.end())
    {
      clusterSize -= 1;
    }

    for (auto atomId : cluster)
    {
      atomClusterTypeVector[atomId] = clusterId;
      atomClusterSizeVector[atomId] = clusterSize;
    }

    clusterSizeVector.emplace_back(clusterSize);

    clusterId++;
  }

  auxiliaryList["clusterType"] = atomClusterTypeVector;
  auxiliaryList["clusterSize"] = atomClusterSizeVector;
}

bool ClusterDynamics::growB2Cluster(size_t atomId,
                                    unordered_set<size_t> &cluster)
{
  if (!visited_.emplace(atomId).second)
    return false;

  bool isB2Center = B2OrderParameter::isB2Ordered(config_, atomId);

  if (!isB2Center)
    return false;

  cluster.insert(atomId);
  // First nearest neighbours
  const auto neighbors = config_.GetNeighborAtomIdVectorOfAtom(atomId, 1);

  for (size_t nId : neighbors)
  {
    cluster.insert(nId); // Include all neighbors regardless of order
    if (B2OrderParameter::isB2Ordered(config_, nId))
    {
      growB2Cluster(nId, cluster); // Recurse on ordered neighbors only
    }
  }
  return true;
}

unordered_set<size_t> ClusterDynamics::mergeIfIntersect(
    const unordered_set<size_t> &set1,
    const unordered_set<size_t> &set2)
{

  for (const size_t &val : set1)
  {
    if (set2.count(val))
    {
      unordered_set<size_t> result = set1;
      result.insert(set2.begin(), set2.end());
      return result;
    }
  }
  return {};
}

void ClusterDynamics::mergeAllClusters()
{
  for (size_t i = 0; i < clusters_.size(); ++i)
  {
    for (size_t j = i + 1; j < clusters_.size();)
    {
      auto combined = mergeIfIntersect(clusters_[i], clusters_[j]);
      if (!combined.empty())
      {
        clusters_.erase(clusters_.begin() + j);
        clusters_.erase(clusters_.begin() + i);
        clusters_.push_back(combined);
        i = static_cast<size_t>(-1);
        break;
      }
      else
      {
        ++j;
      }
    }
  }
}

/*
// Grow a B2 cluster, including all immediate neighbors
bool growB2(const Config &config,
            size_t latticeId,
            unordered_set<size_t> &cluster,
            unordered_set<size_t> &visited)
{
  if (visited.find(latticeId) != visited.end())
  {
    return false;
  }

  visited.emplace(latticeId);

  if (!isB2Ordered(config, latticeId))
  {
    return false;
  }

  cluster.emplace(latticeId);
  cout << "B2-ordered site added: " << latticeId << endl;

  auto neighbors = config.GetNeighborLatticeIdVectorOfLattice(latticeId, 1);
  for (auto nId : neighbors)
  {
    if (cluster.find(nId) == cluster.end())
    {
      cluster.emplace(nId);
      cout << "   Added neighbor: " << nId << endl;
    }
    if (isB2Ordered(config, nId))
    {
      growB2(config, nId, cluster, visited); // Recurse only on B2-ordered neighbors
    }
  }
  return true;
}

// Merge intersecting clusters
unordered_set<size_t> unionIfIntersect(const unordered_set<size_t> &set1,
                                       const unordered_set<size_t> &set2)
{
  for (const size_t &val : set1)
  {
    if (set2.find(val) != set2.end())
    {
      unordered_set<size_t> result = set1;
      result.insert(set2.begin(), set2.end());
      return result;
    }
  }
  return {};
}

// Main clustering function

void clusterGrowth(const Config &config)
{
  size_t numAtoms = config.GetNumAtoms();
  vector<unordered_set<size_t>> allClusters;
  unordered_set<size_t> visited;

  // Step 1: Identify initial clusters
  for (size_t i = 0; i < numAtoms; ++i)
  {
    if (visited.find(i) != visited.end())
    {
      continue;
    }
    unordered_set<size_t> cluster;
    if (growB2(config, i, cluster, visited) && !cluster.empty())
    {
      allClusters.push_back(cluster);
    }
  }

  // Step 2: Merge intersecting clusters
  for (size_t i = 0; i < allClusters.size(); ++i)
  {
    for (size_t j = i + 1; j < allClusters.size();)
    {
      auto combined = unionIfIntersect(allClusters[i], allClusters[j]);
      if (!combined.empty())
      {
        allClusters.erase(allClusters.begin() + j);
        allClusters.erase(allClusters.begin() + i);
        allClusters.push_back(combined);
        i = static_cast<size_t>(-1); // Restart loop
        break;
      }
      else
      {
        ++j;
      }
    }
  }

  // Step 3: Output cluster information
  for (const auto &cluster : allClusters)
  {
    cout << "Cluster size: " << cluster.size() << " : {";
    for (auto id : cluster)
    {
      cout << id << "(" << config.GetElementOfLattice(id) << "), ";
    }
    cout << "}" << endl;
  }

  // Auxiliary list uses atom ids so the properties need to assigned to the atom
  // ids not the lattice ids

  map<string, vector<double>> auxiliary_lists;
  vector<double> clusterSize(numAtoms, -1); // -1 = no cluster
  vector<double> clusterType(numAtoms, -1); // -1 = no cluster

  int clusterTypeId = 1;
  for (const auto &cluster : allClusters)
  {
    for (auto id : cluster)
    {
      // Get the atom id

      size_t atomId = config.GetAtomIdOfLattice(id);
      // clusterSize[id] = static_cast<double>(cluster.size());
      clusterType[atomId] = static_cast<double>(clusterTypeId);
    }
    clusterTypeId++;
  }

  // auxiliary_lists["clusterSize"] = clusterSize;
  auxiliary_lists["clusterType"] = clusterType;

  cout << "Total clusters found: " << allClusters.size() << endl;
  config.GetNeighborAtomIdVectorOfAtom

      Config::WriteConfigExtended("B2ClustersExtendedV4_atom.cfg", config, auxiliary_lists);
}
*/