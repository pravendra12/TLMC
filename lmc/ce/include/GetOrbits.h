#ifndef LMC_CE_INCLUDE_GETORBITS_H_
#define LMC_CE_INCLUDE_GETORBITS_H_

#include <set>
#include <string>
#include <vector>
#include <unordered_map>
#include "Config.h"
#include "ClusterExpansion.h"

using namespace std;

/**
 * @brief Computes the orbits of lattice clusters based on their symmetry and configuration.
 *
 * This function identifies clusters of lattice sites within the given constraints and maps
 * them to their respective orbits based on symmetry. It returns a map where the keys are
 * unique orbit identifiers (as strings) and the values are vectors of clusters (each cluster
 * represented as a vector of site indices).
 *
 * @param config The configuration object containing lattice and simulation parameters.
 * @param maxClusterSize The maximum allowable size of a cluster.
 * @param maxBondOrder The maximum allowable bond order for the clusters.
 * @param equivalentLatticeIdsVector A 2D vector where each sub-vector represents a group of
 *        equivalent sites based on symmetry.
 * @param symmetricallSortedVector A vector of lattice site IDs sorted based on symmetry.
 *
 * @return A map where:
 *         - The key is a string representing a unique orbit identifier (constructed by sorting
 *           and concatenating orbit indices).
 *         - The value is a vector of clusters, where each cluster is represented as a vector
 *           of site indices.
 *
 * @note The function uses helper functions such as `FindClustersWithinAllowedSites` to
 *       identify clusters and assumes the existence of utility functions like `print1DVector`
 *       for debugging purposes.
 *
 * @details
 * - The function first maps lattice site IDs to their corresponding indices and orbits.
 * - It then identifies all clusters within the allowed constraints.
 * - For each cluster, it computes the corresponding orbit vector and generates a unique
 *   key by sorting the orbit indices.
 * - Clusters are grouped into orbits based on their unique keys and stored in the resulting map.
 * - Empty clusters are handled separately and assigned a key of "-1".
 */
map<string, vector<vector<size_t>>> GetOrbits(
    const Config &config,
    const size_t &maxClusterSize,
    const size_t &maxBondOrderOfClusters,
    const vector<size_t> &ssVector,
    const vector<vector<size_t>> &equivalentSitesEncoding);

#endif // LMC_CE_INCLUDE_GETORBITS_H_
