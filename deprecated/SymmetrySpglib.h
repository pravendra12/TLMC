#ifndef LMC_CE_INCLUDE_SYMMETRYSPGLIB_H_
#define LMC_CE_INCLUDE_SYMMETRYSPGLIB_H_

#include "Config.h"
#include <iostream>
#include <vector>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <Eigen/Dense>
#include <spglib.h>
#include "SymmetrySpglib.h"
#include "ClusterExpansion.h"

using namespace std;
using namespace Eigen;

// Hash for Eigen::Matrix3d
struct Matrix3dHash
{
    std::size_t operator()(const Matrix3d &matrix) const
    {
        std::size_t seed = 0;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                boost::hash_combine(seed, matrix(i, j));
        return seed;
    }
};

// Hash for Eigen::Vector3d
struct Vector3dHash2
{
    std::size_t operator()(const Vector3d &vec) const
    {
        std::size_t seed = 0;
        for (int i = 0; i < 3; ++i)
            boost::hash_combine(seed, vec(i));
        return seed;
    }
};

using namespace std;

void get_bcc_symmetry_operations(Config &cfg, vector<size_t> &lce, unordered_set<Matrix3d, Matrix3dHash> &uniqueRotations, unordered_set<Vector3d, Vector3dHash2> &uniqueTranslations, vector<pair<Matrix3d, Vector3d>> &symmetryOperations);

void equivalent_clusters(Config &cfg, vector<size_t> lce);
void equivalent_clusters_triplets(Config &cfg, vector<size_t> lce, size_t latticeId);

vector<pair<Matrix3d, Vector3d>> GetSymmetryOperations(
    const Config &config,
    const unordered_set<size_t> &latticeIdVector,
    const bool debug = false,
    const map<size_t, size_t> &latticeIdToIndexMap = map<size_t, size_t>{},
    const double symprec = 1e-5);

void GetEquivalentClusters(
    const Config &config,
    const unordered_set<size_t> &latticeIdVector,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const double symprec = 1e-5);


vector<size_t> GetCanonicalSortedLatticeSites(
    const Config &config, 
    const pair<size_t, size_t> &latticeIdJumpPair, const size_t &maxBondOrder, const unordered_map<size_t, Eigen::RowVector3d> &referenceLatticeIdHashmap, const vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> &symmetryOperations);

#endif // LMC_CE_INCLUDE_SYMMETRYSPGLIB_H_
