#ifndef LMC_CE_INCLUDE_GETEQUIVALENTCLUSTERS_H_
#define LMC_CE_INCLUDE_GETEQUIVALENTCLUSTERS_H_

#include "Config.h"
#include <iostream>
#include <vector>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <Eigen/Dense>
#include <spglib.h>
#include "SymmetrySpglib.h"
#include "ClusterExpansion.h"
#include "PrintUtility.h"
#include <queue>
#include "UnionFind.h"

using Vec3 = Eigen::Vector3d;

vector<set<vector<size_t>>> GetEquivalentClustersFast(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const double symprec,
    const bool debug,
    const bool isNew = true);

vector<set<vector<size_t>>> GetEquivalentClustersSlow(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet, const double symprec, const bool debug);

vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetOrbitsNew(
  const Config &config, 
  const unordered_set<size_t> &latticeIdSet, 
  const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet, 
  const vector<pair<Matrix3d, Vector3d>> &symmetryOperations, 
  const double symprec = 1e-5, 
  const bool debug = false);

void DebugCluster(const std::vector<size_t> &cluster, const std::unordered_map<size_t, Vec3> &latticeIdToPos,
                  const std::vector<std::pair<Eigen::Matrix3d, Vec3>> &symOps, double scale = 1e8);

void DebugClusterFast(const std::vector<size_t> &cluster, const std::vector<size_t> &latticeIdList,
                      const std::vector<Eigen::Vector3d> &idx_to_frac,
                      const std::vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>> &symOps, double scale = 1e8);

vector<set<vector<size_t>>> GetEquivalentClustersSlowFixed(
    const Config &config, const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const double symprec = 1e-5, const bool debug = true);

vector<set<vector<size_t>>> GetEquivalentClustersFast2(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const vector<pair<Matrix3d, Vector3d>> &symmetryOperations,

    const double symprec, const bool debug, const bool);

vector<set<vector<size_t>>> GetEquivalentClustersByPermutation(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet, const bool debug);

#endif // LMC_CE_INCLUDE_GETEQUIVALENTCLUSTERS_H_