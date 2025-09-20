#ifndef LMC_CE_INCLUDE_CORRELATIONVECTOR_H_
#define LMC_CE_INCLUDE_CORRELATIONVECTOR_H_

#include <Eigen/Dense>
#include "Config.h"
#include "BasisSet.h"
#include "TiledSupercell.h"
#include "LatticeClusterType.hpp"
#include "CorrelationFunction.h"

using namespace std;
using namespace Eigen;

VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters);

// Return correlation vector with a given element assigned to a lattice Id
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const size_t &targetLatticeId,
    const Element &elementToAssign,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters);

VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentOrbitVector);

double GetOrbitCorrelationFunction(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<vector<size_t>> &orbitVector);

double GetOrbitCorrelationFunction(
    const Config &config,
    BasisSet &atomicBasis,
    const size_t &targetLatticeId,
    const Element &elementToAssign,
    const vector<vector<size_t>> &orbitVector);

VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentOrbitVector);

double GetOrbitCorrelationFunction(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<vector<size_t>> &orbitVector);

// TLMC
// CE & LVFE
double GetOrbitCorrelationFunction(
    const TiledSupercell &tiledSupercell,
    const size_t &smallConfigIdx,
    BasisSet &atomicBasis,
    const vector<size_t> &neighborsOfSmallConfig,
    const vector<vector<LatticeSiteEncodedMapping>> &orbitVector);

VectorXd GetCorrelationVector(
    const TiledSupercell &tiledSupercell,
    const LatticeSiteMapping &latticeSite,
    BasisSet &atomicBasis,
    const vector<LatticeSiteEncodedMapping> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters);

double GetOrbitCorrelationFunction(
    const TiledSupercell &tiledSupercell,
    const size_t &smallConfigIdx,
    BasisSet &atomicBasis,
    const LatticeSiteMapping &targetLatticeSite,
    const Element &elementToAssign,
    const vector<size_t> &neighborsOfSmallConfig,
    const vector<vector<LatticeSiteEncodedMapping>> &orbitVector);

VectorXd GetCorrelationVector(
    const TiledSupercell &tiledSupercell,
    const LatticeSiteMapping &latticeSite,
    BasisSet &atomicBasis,
    const LatticeSiteMapping &targetLatticeSite,
    const Element &elementToAssign,
    const vector<LatticeSiteEncodedMapping> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters);

// KRA

VectorXd GetCorrelationVector(
  const TiledSupercell &tiledSupercell, 
  const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair, 
  BasisSet &atomicBasis, const vector<NeighbourOfPair> &canonicalSortedLatticeIds, 
  const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters);

double GetOrbitCorrelationFunction(
    const TiledSupercell &tiledSupercell,
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair,
    BasisSet &atomicBasis,
    const vector<vector<NeighbourOfPair>> &orbitVector);
#endif // LMC_CE_INCLUDE_CORRELATIONVECTOR_H_
