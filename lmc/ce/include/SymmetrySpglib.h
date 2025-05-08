#ifndef LMC_CE_INCLUDE_SYMMETRYSPGLIB_H_
#define LMC_CE_INCLUDE_SYMMETRYSPGLIB_H_

#include "Config.h"
#include "Eigen/Dense"
#include "Constants.hpp"
#include "EncodingUtility.h"
#include "PrintUtility.h"
#include <spglib.h>
#include <Eigen/Dense>
#include <unordered_map>
#include <set>
#include "ClusterExpansion.h"

#include "ClusterExpansion.h"

using namespace std;
using namespace Eigen;


#endif //LMC_CE_INCLUDE_SYMMETRYSPGLIB_H_

std::vector<int> GetSpglibTypes(const Config &config);

void GetSymmetryMatrix(const Config &config);

std::vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>>  GetSymmetryMatrixPure();

LatticeCluster apply_symmetry(const Config &config, const LatticeCluster &latticeCluster);

vector<LatticeCluster> apply_symmetry_cluster(const Config &config, const LatticeCluster &latticeCluster);
