/**************************************************************************************************
 * Copyright (c) 2022-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 6/14/22 12:36 PM                                                                        *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 10/30/23 3:13 PM                                                          *
 **************************************************************************************************/

/*! \file  PotentialEnergyEstimator.cpp
 *  \brief File for the PotentialEnergyEstimator class implementation.
 */

#include "PotentialEnergyEstimator.h"
#include <fstream>
#include <omp.h>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

/*! \brief Convert cluster set to a map with the number of appearance of each cluster type.
 *  \param cluster_type_set : The set of cluster types
 *  \return                 : A map with the number of appearance of each cluster type
 */
static std::unordered_map<ClusterType, size_t, boost::hash<ClusterType>> ConvertSetToHashMap(
    const std::set<ClusterType> &cluster_type_set) {
  std::unordered_map<ClusterType, size_t, boost::hash<ClusterType>> cluster_type_count;
  for (const auto &cluster_type : cluster_type_set) {
    cluster_type_count[cluster_type] = 0;
  }
  return cluster_type_count;
}


/*! \brief Reads the theta values from trained data from json file
 *  \param json_filename Path to the JSON file.
 *  \param json_key Key under which "theta" values are stored.
 *  \return        Eigen::VectorXd containing theta values.
*/
static Eigen::VectorXd ReadVectorFromJson(const std::string &json_filename, const std::string &json_key) {
  std::ifstream ifs(json_filename, std::ifstream::in);
  json all_parameters;
  ifs >> all_parameters;
  Eigen::VectorXd eigen_vector;
  for (const auto &[key, parameters] : all_parameters.items()) {
    if (key == json_key) {
      auto base_theta_json = parameters.at("theta");
      eigen_vector = {};
      for (const auto &theta : base_theta_json) {
        eigen_vector.conservativeResize(eigen_vector.size() + 1);
        eigen_vector(eigen_vector.size() - 1) = theta.get<double>();
      }
    }
  }
  return eigen_vector;
}

PotentialEnergyEstimator::PotentialEnergyEstimator(const std::string &predictor_filename,
                                                   const Config &reference_config,
                                                   const Config &supercell_config,
                                                   const std::set<Element> &element_set,
                                                   size_t max_cluster_size,
                                                   size_t max_bond_order)
    : effective_cluster_interaction_(ReadVectorFromJson(predictor_filename, "Base")),
      element_set_(element_set),
      initialized_cluster_type_set_(
          InitializeClusterTypeSet(reference_config, element_set_, max_cluster_size, max_bond_order)),
          
      lattice_cluster_type_count_(CountLatticeClusterTypes(supercell_config, max_cluster_size, max_bond_order)),
      // all_lattice_hashset_(FindAllLatticeClusters(reference_config, max_cluster_size, max_bond_order, {})),
      max_cluster_size_(max_cluster_size), max_bond_order_(max_bond_order) {

  ///Todo : check the size of the effective_cluster_interaction_ and initialized_cluster_type_set_
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

Eigen::VectorXd PotentialEnergyEstimator::GetEncodeVector(const Config &config) const {
  auto cluster_type_count_hashmap(ConvertSetToHashMap(initialized_cluster_type_set_));
  
  auto all_lattice_hashset = FindAllLatticeClusters(config,3,3,{});

  for (const auto &lattice_cluster : all_lattice_hashset) {
    auto atom_cluster_type = IndentifyAtomClusterType(config, lattice_cluster.GetLatticeIdVector());
    cluster_type_count_hashmap.at(ClusterType(atom_cluster_type, lattice_cluster.GetClusterType()))++;
  }
  Eigen::VectorXd encode_vector(initialized_cluster_type_set_.size());
  int idx = 0;
  for (const auto &cluster_type : initialized_cluster_type_set_) {
    // Count of Cluster Types for given configuration
    auto count_bond = static_cast<double>(cluster_type_count_hashmap.at(cluster_type));
    // Count of Cluster Types for normalization
    auto total_bond = static_cast<double>(lattice_cluster_type_count_.at(cluster_type.lattice_cluster_type_));

    encode_vector(idx) = count_bond/total_bond;

    std::cout << cluster_type << " : " << count_bond << " : " << total_bond << std::endl;
    ++idx;

  }

  return encode_vector;
}


Eigen::VectorXd PotentialEnergyEstimator::GetEncodeVectorOfCluster(const Config &config, std::vector<size_t> cluster) const {
  auto cluster_type_count_hashmap(ConvertSetToHashMap(initialized_cluster_type_set_));

  auto all_lattice_clusters  = FindAllLatticeClusters(config, max_cluster_size_, max_bond_order_, cluster);

  for (auto &lattice_cluster : all_lattice_clusters) {
    auto atom_cluster_type = IndentifyAtomClusterType(config, lattice_cluster.GetLatticeIdVector());
    cluster_type_count_hashmap.at(ClusterType(atom_cluster_type, lattice_cluster.GetClusterType()))++;
  }

  Eigen::VectorXd encode_vector_cluster(initialized_cluster_type_set_.size());
  int idx = 0;
  for (const auto &cluster_type : initialized_cluster_type_set_) {
    // Count of Cluster Types for given configuration
    auto count_bond = static_cast<double>(cluster_type_count_hashmap.at(cluster_type));
    // Count of Cluster Types for normalization
    auto total_bond = static_cast<double>(lattice_cluster_type_count_.at(cluster_type.lattice_cluster_type_));

    encode_vector_cluster(idx) = count_bond/total_bond;

    //std::cout << cluster_type << " : " << count_bond << " : " << total_bond << std::endl;
    ++idx;

  }

  return encode_vector_cluster;
}

double PotentialEnergyEstimator::GetEnergy(const Config &config) const {
    auto encode_vector = GetEncodeVector(config);

    // std::cout << "I am here ..... " << std::endl;

    // Eigen::VectorXd energy_vector = effective_cluster_interaction_.array()*encode_vector.array();

    Eigen::VectorXd energy_vector(encode_vector.size());
    for (int i = 0; i < encode_vector.size(); ++i) {
        energy_vector[i] = effective_cluster_interaction_[i] * encode_vector[i];
    }

    return energy_vector.sum();
}

double PotentialEnergyEstimator::GetEnergyOfCluster(const Config &config, const std::vector<size_t> &cluster) const {
  auto encode_vector_cluster = GetEncodeVectorOfCluster(config, cluster);
  
  // Eigen::VectorXd energy_vector = effective_cluster_interaction_.array()*encode_vector_cluster.array();

  Eigen::VectorXd energy_vector(encode_vector_cluster.size());
  for (int i = 0; i < encode_vector_cluster.size(); ++i) {
    energy_vector[i] = effective_cluster_interaction_[i] * encode_vector_cluster[i];
  }
   
  return energy_vector.sum();

}

double PotentialEnergyEstimator::GetDe(Config &config, const std::pair<size_t, size_t> &lattice_id_pair) const {
  
  if (config.GetElementOfLattice(lattice_id_pair.first) == config.GetElementOfLattice(lattice_id_pair.second)) {
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