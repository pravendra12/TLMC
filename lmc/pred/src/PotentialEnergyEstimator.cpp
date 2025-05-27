/**************************************************************************************************
 * Copyright (c) 2022-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 6/14/22 12:36 PM                                                                        *
 * @Last Modified by: pravendra12                                                                  *
 * @Last Modified time: 10/30/23 3:13 PM                                                          *
 **************************************************************************************************/

/*! \file  PotentialEnergyEstimator.cpp
 *  \brief File for the PotentialEnergyEstimator class implementation.
 */

#include "PotentialEnergyEstimator.h"

/*! \brief Convert cluster set to a map with the number of appearance of each cluster type.
 *  \param cluster_type_set : The set of cluster types
 *  \return                 : A map with the number of appearance of each cluster type
 */
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

PotentialEnergyEstimator::PotentialEnergyEstimator(
    const std::string &predictor_filename,
    const Config &reference_config,
    const Config &supercell_config,
    const std::set<Element> &element_set) : max_cluster_size_(
                                                ReadParameterFromJson(
                                                predictor_filename,
                                                "maxClusterSizeCE")),
                                            max_bond_order_(
                                                ReadParameterFromJson(
                                                    predictor_filename,
                                                    "maxBondOrderCE")),
                                            beta_ce_(
                                              ReadParametersFromJson(
                                                predictor_filename, 
                                                "ce", "beta_ce")),
                                            element_set_(element_set),
                                            initialized_cluster_type_set_(
                                                InitializeClusterTypeSet(reference_config,
                                                                         element_set_,
                                                                         max_cluster_size_,
                                                                         max_bond_order_)),
                                            lattice_cluster_type_count_(
                                                CountLatticeClusterTypes(supercell_config,
                                                                         max_cluster_size_,
                                                                         max_bond_order_))

{

  cout << "Max Bond Order for CE: " << max_bond_order_ << endl;
  cout << "Max Cluster Size for CE: " << max_cluster_size_ << endl;


  /// Todo : check the size of the effective_cluster_interaction_ and initialized_cluster_type_set_
  // if (initialized_cluster_type_set_.size() != static_cast<size_t>(effective_cluster_interaction_.size())) {
  //   // here 1 is for void cluster
  //   throw std::invalid_argument(
  //       "The size of the ECI vector is not compatible with the number of types of clusters. They are "
  //           + std::to_string(effective_cluster_interaction_.size()) + " and "
  //           + std::to_string(initialized_cluster_type_set_.size()) + " respectively.");
  // }
}
PotentialEnergyEstimator::~PotentialEnergyEstimator() = default;

// This is the function used to generate the encode_vector which has the count
// corresponding to each cluster

Eigen::VectorXd PotentialEnergyEstimator::GetEncodeVector(const Config &config) const
{
  auto cluster_type_count_hashmap(ConvertSetToHashMap(initialized_cluster_type_set_));

  auto all_lattice_hashset = FindAllLatticeClusters(config, max_cluster_size_, max_bond_order_, {});

  for (const auto &lattice_cluster : all_lattice_hashset)
  {
    auto atom_cluster_type = IdentifyAtomClusterType(config, lattice_cluster.GetLatticeIdVector());
    cluster_type_count_hashmap.at(ClusterType(atom_cluster_type, lattice_cluster.GetClusterType()))++;
  }
  Eigen::VectorXd encode_vector(initialized_cluster_type_set_.size());
  int idx = 0;
  for (const auto &cluster_type : initialized_cluster_type_set_)
  {
    // Count of Cluster Types for given configuration
    auto count_bond = static_cast<double>(cluster_type_count_hashmap.at(cluster_type));
    // Count of Cluster Types for normalization
    auto total_bond = static_cast<double>(lattice_cluster_type_count_.at(cluster_type.lattice_cluster_type_));

    encode_vector(idx) = count_bond / total_bond;

    std::cout << cluster_type << " : " << count_bond << " : " << total_bond << std::endl;
    ++idx;
  }

  return encode_vector;
}

Eigen::VectorXd
PotentialEnergyEstimator::GetEncodeVectorOfCluster(const Config &config,
                                                   const std::vector<size_t> &cluster) const
{
  auto cluster_type_count_hashmap(ConvertSetToHashMap(initialized_cluster_type_set_));

  auto all_lattice_clusters = FindAllLatticeClusters(config, max_cluster_size_, max_bond_order_, cluster);

  for (auto &lattice_cluster : all_lattice_clusters)
  {
    auto atom_cluster_type = IdentifyAtomClusterType(config, lattice_cluster.GetLatticeIdVector());
    cluster_type_count_hashmap.at(ClusterType(atom_cluster_type, lattice_cluster.GetClusterType()))++;
  }

  Eigen::VectorXd encode_vector_cluster(initialized_cluster_type_set_.size());
  int idx = 0;
  for (const auto &cluster_type : initialized_cluster_type_set_)
  {
    // Count of Cluster Types for given configuration
    auto count_bond = static_cast<double>(cluster_type_count_hashmap.at(cluster_type));
    // Count of Cluster Types for normalization
    auto total_bond = static_cast<double>(lattice_cluster_type_count_.at(cluster_type.lattice_cluster_type_));

    encode_vector_cluster(idx) = count_bond / total_bond;

    //std::cout << cluster_type << " : " << count_bond << " : " << total_bond << std::endl;

    ++idx;
  }

  return encode_vector_cluster;
}

double
PotentialEnergyEstimator::GetEnergy(const Config &config) const
{
  auto encode_vector = GetEncodeVector(config);

  // double energy = adjusted_beta_ce_.dot(encode_vector) + adjusted_intercept_ce_;

  double energy = beta_ce_.dot(encode_vector);

  return energy;
}

double
PotentialEnergyEstimator::GetEnergyOfCluster(const Config &config,
                                             const std::vector<size_t> &cluster) const
{
  Eigen::VectorXd encode_vector_cluster = GetEncodeVectorOfCluster(config,
                                                                   cluster);

  // double energy_cluster = adjusted_beta_ce_.dot(encode_vector_cluster) +
  //                         adjusted_intercept_ce_;
  double energy_cluster = beta_ce_.dot(encode_vector_cluster);
  return energy_cluster;
}

double
PotentialEnergyEstimator::GetDe(Config &config,
                                const std::pair<size_t, size_t> &lattice_id_pair) const
{
  if (config.GetElementOfLattice(lattice_id_pair.first) == config.GetElementOfLattice(lattice_id_pair.second))
  {
    return 0;
  }

  // Energy Before Swap
  auto E_before_swap = GetEnergyOfCluster(config, {lattice_id_pair.first, lattice_id_pair.second});

  // Swapping the Elements
  config.LatticeJump(lattice_id_pair);

  // Energy After Swap
  auto E_after_swap = GetEnergyOfCluster(config, {lattice_id_pair.first, lattice_id_pair.second});

  // Going back to Original Config
  config.LatticeJump(lattice_id_pair);

  auto dE = E_after_swap - E_before_swap;

  return dE;
}

double PotentialEnergyEstimator::GetDeThreadSafe(
    const Config &config,
    const std::pair<size_t, size_t> &lattice_id_pair) const
{
  const size_t id1 = lattice_id_pair.first;  // Vacancy (X)
  const size_t id2 = lattice_id_pair.second; // Migrating atom

  if (config.GetElementOfLattice(id1) == config.GetElementOfLattice(id2))
  {
    // cout << "0" << endl;
    return 0;
  }

  // Step 1: Find all clusters affected by the jump pair
  std::vector<size_t> cluster = {id1, id2};

  auto start = chrono::high_resolution_clock::now();
  auto all_lattice_clusters = FindAllLatticeClusters(config, max_cluster_size_, max_bond_order_, cluster);
  auto end = chrono::high_resolution_clock::now();
  auto duration = duration_cast<chrono::microseconds>(end - start);
  // cout << "Time to compute FindAllLatticeClusters: " << duration.count() << " microseconds" << endl;

  // Step 2: Compute encoding contributions before and after the swap
  auto cluster_type_count_hashmap_before = ConvertSetToHashMap(initialized_cluster_type_set_);
  auto cluster_type_count_hashmap_after = cluster_type_count_hashmap_before;

  for (const auto &lattice_cluster : all_lattice_clusters)
  {
    auto lattice_ids = lattice_cluster.GetLatticeIdVector();

    // Before swap: Use current config state
    auto atom_cluster_type_before = IdentifyAtomClusterType(config, lattice_ids);
    auto cluster_type_before = ClusterType(atom_cluster_type_before, lattice_cluster.GetClusterType());
    cluster_type_count_hashmap_before.at(cluster_type_before)++;

    // After swap: Simulate the swap locally
    std::vector<Element> swapped_elements;
    for (auto id : lattice_ids)
    {
      if (id == id1)
      {
        swapped_elements.push_back(config.GetElementOfLattice(id2)); // id1 gets id2's element
      }
      else if (id == id2)
      {
        // swapped_elements.push_back(Element("X")); // id2 gets id1's element
        swapped_elements.push_back(config.GetElementOfLattice(id1));
      }
      else
      {
        swapped_elements.push_back(config.GetElementOfLattice(id));
      }
    }
    auto atom_cluster_type_after = AtomClusterType(swapped_elements);
    auto cluster_type_after = ClusterType(atom_cluster_type_after, lattice_cluster.GetClusterType());
    cluster_type_count_hashmap_after.at(cluster_type_after)++;
  }

  // Step 3: Convert to encoding vectors
  Eigen::VectorXd encode_before(initialized_cluster_type_set_.size());
  Eigen::VectorXd encode_after(initialized_cluster_type_set_.size());
  int idx = 0;
  for (const auto &cluster_type : initialized_cluster_type_set_)
  {
    auto count_before = static_cast<double>(cluster_type_count_hashmap_before.at(cluster_type));
    auto count_after = static_cast<double>(cluster_type_count_hashmap_after.at(cluster_type));
    auto total_bond = static_cast<double>(lattice_cluster_type_count_.at(cluster_type.lattice_cluster_type_));
    encode_before(idx) = count_before / total_bond;
    encode_after(idx) = count_after / total_bond;
    ++idx;

    // cout << cluster_type << " : " << count_before << " : " << count_after << endl;

  }

  cout <<  "encodeBefore : " << encode_before.transpose() << endl;
  cout <<  "encodeAfter : " << encode_after.transpose() << endl;


  // Step 4: Compute energy difference

  // double E_before_swap = beta_ce_.dot(encode_before);
  // double E_after_swap = beta_ce_.dot(encode_after);
  //
  // double dE = E_after_swap - E_before_swap;

  VectorXd encode_diff = encode_after - encode_before;
  double dE = beta_ce_.dot(encode_diff);

  // cout << dE << endl;

  return dE;
}