/*******************************************************************************
 * Copyright (c) 2022-2023. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 6/14/22 12:36 PM
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! @file  PotentialEnergyEstimator.cpp
 *  @brief File for the PotentialEnergyEstimator class implementation.
 */

#include "PotentialEnergyEstimator.h"

/*! @brief Convert cluster set to a map with the number of appearance of each cluster type.
 *  @param clusterTypeSet  The set of cluster types
 *  @return                A map with the number of appearance of each cluster type
 */
static unordered_map<ClusterType, size_t, boost::hash<ClusterType>> ConvertSetToHashMap(
    const set<ClusterType> &clusterTypeSet)
{
  unordered_map<ClusterType, size_t, boost::hash<ClusterType>> clusterTypeCount;
  for (const auto &clusterType : clusterTypeSet)
  {
    clusterTypeCount[clusterType] = 0;
  }
  return clusterTypeCount;
}

PotentialEnergyEstimator::PotentialEnergyEstimator(
    const string &predictorFilename,
    const Config &referenceConfig,
    const Config &supercellConfig,
    const set<Element> &element_set) : maxClusterSize_(ReadParameterFromJson(predictorFilename,
                                                                             "maxClusterSizeCE")),
                                       maxBondOrder_(
                                           ReadParameterFromJson(
                                               predictorFilename,
                                               "maxBondOrderCE")),
                                       betaCE_(
                                           ReadParametersFromJson(
                                               predictorFilename,
                                               "ce", "beta_ce")),
                                       elementSet_(element_set),

                                       initializedClusterTypeSet_(
                                           InitializeClusterTypeSet(
                                               referenceConfig,
                                               elementSet_,
                                               maxClusterSize_,
                                               maxBondOrder_)),
                                       clusterTypeCountHashMap_(
                                           ConvertSetToHashMap(
                                               initializedClusterTypeSet_)),
                                       latticeClusterTypeCount_(
                                           CountLatticeClusterTypes(
                                               supercellConfig,
                                               maxClusterSize_,
                                               maxBondOrder_))

{

  cout << "Max Bond Order for CE: " << maxBondOrder_ << endl;
  cout << "Max Cluster Size for CE: " << maxClusterSize_ << endl;

  if (initializedClusterTypeSet_.size() != static_cast<size_t>(betaCE_.size()))
  {
    throw std::invalid_argument(
        "Error in PotentialEnergyEstimator: The size of 'betaCE_' (" + std::to_string(betaCE_.size()) + ") is not compatible with the number of initialized cluster types in 'initializedClusterTypeSet_' (" + std::to_string(initializedClusterTypeSet_.size()) + "). Please ensure they are consistent.");
  }
}

PotentialEnergyEstimator::~PotentialEnergyEstimator() = default;

Eigen::VectorXd PotentialEnergyEstimator::GetEncodeVector(const Config &config) const
{
  auto clusterTypeCountHashMap = clusterTypeCountHashMap_;

  auto allLatticeClusterSet = FindAllLatticeClusters(config, maxClusterSize_, maxBondOrder_, {});

  for (const auto &latticeCluster : allLatticeClusterSet)
  {
    auto atomClusterType = IdentifyAtomClusterType(config, latticeCluster.GetLatticeIdVector());
    clusterTypeCountHashMap.at(ClusterType(atomClusterType, latticeCluster.GetClusterType()))++;
  }
  Eigen::VectorXd encodeVector(initializedClusterTypeSet_.size());
  int idx = 0;
  for (const auto &clusterType : initializedClusterTypeSet_)
  {
    // Count of Cluster Types for given configuration
    auto count_bond = static_cast<double>(clusterTypeCountHashMap.at(clusterType));
    // Count of Cluster Types for normalization
    auto total_bond = static_cast<double>(latticeClusterTypeCount_.at(clusterType.lattice_cluster_type_));

    encodeVector(idx) = count_bond / total_bond;
    ++idx;
  }

  return encodeVector;
}

Eigen::VectorXd PotentialEnergyEstimator::GetEncodeVectorOfCluster(
    const Config &config,
    const vector<size_t> &cluster) const
{
  auto clusterTypeCountHashMap = clusterTypeCountHashMap_;

  auto allLatticeClusterSet = FindAllLatticeClusters(config, maxClusterSize_, maxBondOrder_, cluster);

  for (auto &latticeCluster : allLatticeClusterSet)
  {
    auto atomClusterType = IdentifyAtomClusterType(config, latticeCluster.GetLatticeIdVector());
    clusterTypeCountHashMap.at(ClusterType(atomClusterType, latticeCluster.GetClusterType()))++;
  }

  Eigen::VectorXd encodeVectorCluster(initializedClusterTypeSet_.size());
  int idx = 0;
  for (const auto &clusterType : initializedClusterTypeSet_)
  {
    // Count of Cluster Types for given configuration
    auto count_bond = static_cast<double>(clusterTypeCountHashMap.at(clusterType));
    // Count of Cluster Types for normalization
    auto total_bond = static_cast<double>(latticeClusterTypeCount_.at(clusterType.lattice_cluster_type_));

    encodeVectorCluster(idx) = count_bond / total_bond;
    ++idx;
  }

  return encodeVectorCluster;
}

double PotentialEnergyEstimator::GetEnergy(const Config &config) const
{
  auto encodeVector = GetEncodeVector(config);

  double energy = betaCE_.dot(encodeVector);

  return energy;
}

// Need to change the name of this
double PotentialEnergyEstimator::GetEnergyOfCluster(const Config &config,
                                                    const vector<size_t> &cluster) const
{
  Eigen::VectorXd encodeVectorCluster = GetEncodeVectorOfCluster(config,
                                                                 cluster);

  double energyOfCluster = betaCE_.dot(encodeVectorCluster);
  return energyOfCluster;
}

Eigen::VectorXd PotentialEnergyEstimator::GetEncodeVectorWithinAllowedSites(
    const Config &config,
    const vector<size_t> &allowedSites) const
{
  auto clusterTypeCountHashMap = clusterTypeCountHashMap_;

  auto allLatticeClusterSet = FindClustersWithinAllowedSites(config, maxClusterSize_, maxBondOrder_, allowedSites);
  for (auto &latticeCluster : allLatticeClusterSet)
  {
    auto atomClusterType = IdentifyAtomClusterType(config, latticeCluster.GetLatticeIdVector());
    clusterTypeCountHashMap.at(ClusterType(atomClusterType, latticeCluster.GetClusterType()))++;
  }

  Eigen::VectorXd encodeVectorCluster(initializedClusterTypeSet_.size());
  int idx = 0;
  for (const auto &clusterType : initializedClusterTypeSet_)
  {
    // Count of Cluster Types for given configuration
    auto count_bond = static_cast<double>(clusterTypeCountHashMap.at(clusterType));
    // Count of Cluster Types for normalization
    auto total_bond = static_cast<double>(latticeClusterTypeCount_.at(clusterType.lattice_cluster_type_));

    encodeVectorCluster(idx) = count_bond / total_bond;

    // cout << clusterType.lattice_cluster_type_ << " " << clusterType.atom_cluster_type_ << " : " << count_bond << endl;

    ++idx;
  }

  return encodeVectorCluster;
}

double PotentialEnergyEstimator::GetEnergyOfClusterWithinAllowedSites(
    const Config &config,
    const vector<size_t> &allowedSites) const
{
  // allowed sites would contain the lattice Ids
  Eigen::VectorXd encodeVectorCluster = GetEncodeVectorWithinAllowedSites(config,
                                                                          allowedSites);

  double energyOfCluster = betaCE_.dot(encodeVectorCluster);
  return energyOfCluster;
}

double PotentialEnergyEstimator::GetDeSwap(Config &config,
                                           const pair<size_t, size_t> &latticeIdPair) const
{
  if (config.GetElementOfLattice(latticeIdPair.first) == config.GetElementOfLattice(latticeIdPair.second))
  {
    return 0;
  }

  // Energy Before Swap
  auto energyBeforeSwap = GetEnergyOfCluster(config, {latticeIdPair.first, latticeIdPair.second});

  // Swapping the Elements
  config.LatticeJump(latticeIdPair);

  // Energy After Swap
  auto energyAfterSwap = GetEnergyOfCluster(config, {latticeIdPair.first, latticeIdPair.second});

  // Going back to Original Config
  config.LatticeJump(latticeIdPair);

  auto dE = energyAfterSwap - energyBeforeSwap;

  return dE;
}

double PotentialEnergyEstimator::GetDeMigration(
    const Config &config,
    const pair<size_t, size_t> &latticeIdPair) const
{
  const size_t id1 = latticeIdPair.first;  // Vacancy (X)
  const size_t id2 = latticeIdPair.second; // Migrating atom

  if (config.GetElementOfLattice(id1) == config.GetElementOfLattice(id2))
  {
    return 0;
  }

  // Find all clusters affected by the jump pair
  vector<size_t> cluster = {id1, id2};

  auto allLatticeClusterSet = FindAllLatticeClusters(config,
                                                     maxClusterSize_,
                                                     maxBondOrder_,
                                                     cluster);

  // Compute encoding contributions before and after the swap
  auto clusterTypeCountHashMapBefore = clusterTypeCountHashMap_;
  auto clusterTypeCountHashMapAfter = clusterTypeCountHashMap_;

  for (const auto &latticeCluster : allLatticeClusterSet)
  {
    auto lattice_ids = latticeCluster.GetLatticeIdVector();

    // Before swap: Use current config state
    auto atomClusterTypeBefore = IdentifyAtomClusterType(config, lattice_ids);
    auto clusterTypeBefore = ClusterType(atomClusterTypeBefore, latticeCluster.GetClusterType());
    clusterTypeCountHashMapBefore.at(clusterTypeBefore)++;

    // After swap: Simulate the swap locally
    vector<Element> swappedElements;
    for (auto id : lattice_ids)
    {
      if (id == id1)
      {
        // id1 gets id2's element
        swappedElements.push_back(config.GetElementOfLattice(id2));
      }
      else if (id == id2)
      {
        swappedElements.push_back(config.GetElementOfLattice(id1));
      }
      else
      {
        swappedElements.push_back(config.GetElementOfLattice(id));
      }
    }
    auto atomClusterTypeAfter = AtomClusterType(swappedElements);
    auto clusterTypeAfter = ClusterType(atomClusterTypeAfter, latticeCluster.GetClusterType());
    clusterTypeCountHashMapAfter.at(clusterTypeAfter)++;
  }

  // Convert to encoding vectors
  Eigen::VectorXd encodeBefore(initializedClusterTypeSet_.size());
  Eigen::VectorXd encodeAfter(initializedClusterTypeSet_.size());

  int idx = 0;
  for (const auto &clusterType : initializedClusterTypeSet_)
  {
    auto count_before = static_cast<double>(clusterTypeCountHashMapBefore.at(clusterType));
    auto count_after = static_cast<double>(clusterTypeCountHashMapAfter.at(clusterType));
    auto total_bond = static_cast<double>(latticeClusterTypeCount_.at(clusterType.lattice_cluster_type_));

    encodeBefore(idx) = count_before / total_bond;
    encodeAfter(idx) = count_after / total_bond;
    ++idx;
  }

  // Compute energy difference

  VectorXd encode_diff = encodeAfter - encodeBefore;
  double dE = betaCE_.dot(encode_diff);

  return dE;
}