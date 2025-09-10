#ifndef LMC_CE_INCLUDE_SYMMETRY_H_
#define LMC_CE_INCLUDE_SYMMETRY_H_

#include "Config.h"
#include "UnionFind.h"
#include "Constants.hpp"
#include "PrintUtility.h"

#include <cmath>
#include <array>
#include <utility>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace std;
using namespace Eigen;

struct Vector3dHash
{
  std::size_t operator()(const Vector3d &v) const
  {
    return std::hash<double>()(v.x()) ^ std::hash<double>()(v.y()) ^ std::hash<double>()(v.z());
  }
};

vector<vector<size_t>> GetEquivalentSitesUnder3BarSymmetry(
    const Config &config,
    const pair<size_t, size_t> latticeIdJumpPair,
    const size_t maxBondOrder);

vector<size_t> GetSSVector3FSymmetry(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair,
    const size_t maxBondOrder);

vector<vector<size_t>> GetEquivalentSiteEncoding3BarSymmetry(
  const Config &config, 
  size_t maxBondOrder);

#endif // LMC_CE_INCLUDE_SYMMETRY_H_
