#ifndef LMC_CE_INCLUDE_GETEQUIVALENTCLUSTERS_H_
#define LMC_CE_INCLUDE_GETEQUIVALENTCLUSTERS_H_

#include <queue>
#include <vector>
#include <spglib.h>
#include <iostream>
#include <Eigen/Dense>
#include <unordered_set>
#include <boost/functional/hash.hpp>

#include "Config.h"
#include "UnionFind.h"
#include "Constants.hpp"
#include "PrintUtility.h"
#include "SymmetrySpglib.h"
#include "ClusterExpansion.h"

using namespace Eigen;
using namespace std;

// Helper: small struct to hash vector<int64_t>
struct VecInt64Hash
{
  size_t operator()(const vector<int64_t> &v) const noexcept
  {
    size_t seed = v.size();
    for (auto x : v)
    {
      // combine, similar to boost::hash_combine
      seed ^= std::hash<int64_t>{}(x) + 0x9e3779b97f4a7c15ULL // magic constant (Knuth's golden ratio)
              + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};


// Returns clusters grouped based on symmetry
// Symmetry Operation are determined by the latticeIdSet
// For eg. for a single site and its neighbours
//         Symmetry operation will be m-3m
//
//         Whereas for the a local environment of first nn latticeIdPair
//         Symmetry operation will be -3m

vector<set<vector<size_t>>> GetEquivalentClusters(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const bool debug = false,
    const double symprec = constants::SYMPREC);

// Returns the encoded orbits given the sorted lattice Ids and group equivalentClusters
vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetEncodedOrbits(
  const Config &config, 
  const vector<size_t> &sortedNNLatticeIdVector, 
  const vector<set<vector<size_t>>> &equivalentClusters, 
  const bool debug = false);

/**
 * @brief Returns the Local Encoded Orbits for a lattice site which defines the local 
 *        environment around that site. Following steps are used for the same:
 *        1. Get all lattice clusters center around the latticeId.
 *        2. Group these cluster based on the symmetry operation of its neighours.
 *        3. Remove the latticeId from all the cluster, leaving only those clusters
 *           which uniquely define the local surrounding around a site
 *        4. Finally the reduced clusters are encoded using sorted neighbour latticeIds
 * 
 * NOTE : To be consistent, since the site is being removed hence newClusterSize = maxClusterSize + 1
 * 
 * @param config 
 * @param maxBondOrder 
 * @param maxClusterSize 
 * @param debug 
 * @param symprec 
 * @return vector<pair<vector<vector<size_t>>, LatticeClusterType>> 
 */
vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetLocalEncodedOrbitsForSite(
  const Config &config, 
  const size_t &maxBondOrder, 
  const size_t &maxClusterSize, 
  const bool debug = false, 
  const double symprec = constants::SYMPREC);

  /**
 * @brief Returns the Local Encoded Orbits for around a pair which defines the local 
 *        environment around it. Following steps are used for the same:
 *        1. Get all lattice clusters center around the latticeIds in the jumpPair.
 *        2. Group these cluster based on the symmetry operation of its neighours.
 *        3. Remove the latticeId from all the cluster, leaving only those clusters
 *           which uniquely define the local surrounding around a site
 *        4. Finally the reduced clusters are encoded using sorted neighbour latticeIds
 * 
 * NOTE : To be consistent, since the site is being removed hence newClusterSize = maxClusterSize + 1
 * 
 * @param config 
 * @param maxBondOrder 
 * @param maxClusterSize 
 * @param referenceJumpDirection Direction along which the neighbouring ids will be sorted.
 * @param debug 
 * @param symprec 
 * @return vector<pair<vector<vector<size_t>>, LatticeClusterType>> 
 */
vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetLocalEncodedOrbitsForPair(
  const Config &config, 
  const size_t &maxBondOrder, 
  const size_t &maxClusterSize, 
  const Vector3d &referenceJumpDirection, 
  const bool debug = false, 
  const double symprec = constants::SYMPREC);

#endif // LMC_CE_INCLUDE_GETEQUIVALENTCLUSTERS_H_