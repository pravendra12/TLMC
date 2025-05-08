/**
 * @file ClusterExpansion.h
 * @brief Provides declarations for ClusterExpansion methods used in lattice Monte Carlo simulations.
 *
 * This header defines functions for identifying cluster types, initializing cluster type sets,
 * and finding clusters in lattice configurations. It supports operations on lattice clusters,
 * atom clusters, and general clusters, enabling efficient representation and manipulation
 * of atomic arrangements in simulations.
 *
 * @version 0.1
 * @date 2025-05-04
 * @author Zhucong Xi
 * @last_modified_by pravendra
 * @last_modified_time 2024-11-26
 * @copyright Copyright (c) 2023-2025
 */

#ifndef LMC_CE_INCLUDE_CLUSTEREXPANSION_H_
#define LMC_CE_INCLUDE_CLUSTEREXPANSION_H_

#include <set>
#include <numeric>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include "Eigen/Dense"
#include "LatticeCluster.hpp"
#include "AtomCluster.hpp"
#include "ClusterType.hpp"

/**
 * @brief Identifies the type of a lattice cluster.
 *
 * @param reference_config The configuration to analyze.
 * @param cluster The input cluster, represented by lattice IDs.
 * @return The lattice type of the cluster.
 */
LatticeClusterType IdentifyLatticeClusterType(const Config &reference_config,
                                               const std::vector<size_t> &cluster);

/**
 * @brief Identifies the type of an atom cluster.
 *
 * @param reference_config The configuration to analyze.
 * @param cluster The input cluster, represented by lattice IDs.
 * @return The atom type of the cluster.
 */
AtomClusterType IdentifyAtomClusterType(const Config &reference_config,
                                         const std::vector<size_t> &cluster);

/**
 * @brief Identifies the type of a general cluster.
 *
 * @param reference_config The configuration to analyze.
 * @param cluster The input cluster.
 * @return The type of the cluster.
 */
ClusterType IdentifyClusterType(const Config &reference_config,
                                 const std::vector<size_t> &cluster);

/**
 * @brief Initializes a set of all available lattice cluster types.
 *
 * @param reference_config The configuration to analyze.
 * @param max_cluster_size The maximum size of clusters.
 * @param max_bond_order The cutoff bond distance order applied.
 * @return A set of lattice cluster types.
 */
std::set<LatticeClusterType> InitializeLatticeClusterTypeSet(const Config &reference_config,
                                                             size_t max_cluster_size,
                                                             size_t max_bond_order);

/**
 * @brief Initializes a set of all available atom cluster types.
 *
 * @param element_set The set of elements.
 * @param max_cluster_size The maximum size of clusters.
 * @return A set of atom cluster types.
 */
std::set<AtomClusterType> InitializeAtomClusterTypeSet(const std::set<Element> &element_set,
                                                       size_t max_cluster_size);

/**
 * @brief Initializes a set of all available cluster types.
 *
 * @param reference_config The configuration to analyze.
 * @param element_set The set of elements.
 * @param max_cluster_size The maximum size of clusters.
 * @param max_bond_order The cutoff bond distance order applied.
 * @return A set of cluster types.
 */
std::set<ClusterType> InitializeClusterTypeSet(const Config &reference_config,
                                               const std::set<Element> &element_set,
                                               size_t max_cluster_size,
                                               size_t max_bond_order);

/**
 * @brief Finds all lattice clusters related to a given lattice ID vector.
 *
 * @param reference_config The configuration to analyze.
 * @param max_cluster_size The maximum size of clusters.
 * @param max_bond_order The cutoff bond distance order applied.
 * @param lattice_id_vector The vector of lattice IDs.
 * @return A hashset of lattice clusters.
 */
std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>> FindAllLatticeClusters(
    const Config &reference_config,
    size_t max_cluster_size,
    size_t max_bond_order,
    const std::vector<size_t> &lattice_id_vector);

/**
 * @brief Counts the occurrences of each lattice cluster type in the configuration.
 *
 * @param reference_config The configuration to analyze.
 * @param max_cluster_size The maximum size of clusters.
 * @param max_bond_order The cutoff bond distance order applied.
 * @return A hashmap mapping lattice cluster types to their counts.
 */
std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>>
CountLatticeClusterTypes(const Config &reference_config,
                         const size_t max_cluster_size,
                         const size_t max_bond_order);

/**
 * @brief Counts the occurrences of lattice sites in the configuration.
 *
 * @param reference_config The configuration to analyze.
 * @param max_cluster_size The maximum size of clusters.
 * @param max_bond_order The cutoff bond distance order applied.
 * @return A hashmap mapping lattice cluster types to their counts.
 */
std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>>
CountLatticeSite(const Config &reference_config,
                 const size_t max_cluster_size,
                 const size_t max_bond_order);

/**
 * @brief Creates a hashmap of lattice sites and their associated clusters.
 *
 * @param reference_config The configuration to analyze.
 * @param max_cluster_size The maximum size of clusters.
 * @param max_bond_order The cutoff bond distance order applied.
 * @return A hashmap mapping lattice cluster types to their associated clusters.
 */
std::unordered_map<LatticeClusterType,
                   std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>>,
                   boost::hash<LatticeClusterType>>
LatticeSiteHashMap(const Config &reference_config,
                   const size_t max_cluster_size,
                   const size_t max_bond_order);

/**
 * @brief Finds all clusters that can be formed using allowed lattice sites.
 *
 * Constructs clusters from a reference configuration using only the specified lattice sites.
 *
 * @param config The reference lattice configuration.
 * @param maxClusterSize The maximum number of sites allowed in each cluster.
 * @param maxBondOrder The maximum neighbor shell (bond distance) used when forming clusters.
 * @param allowedLatticeSites A list of lattice site IDs allowed in the clusters.
 * @return A set of clusters representing the local chemical environment.
 */
std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>> FindClustersWithinAllowedSites(
    const Config &config,
    const size_t maxClusterSize,
    const size_t maxBondOrder,
    const std::vector<size_t> &allowedLatticeSites);

#endif // LMC_CE_INCLUDE_CLUSTEREXPANSION_H_
