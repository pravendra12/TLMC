#ifndef LMC_CE_INCLUDE_CORRELATIONVECTOR_H_
#define LMC_CE_INCLUDE_CORRELATIONVECTOR_H_


#include <Eigen/Dense>
#include "Config.h"
#include "BasisSet.h"
#include "LatticeClusterType.hpp"
#include "CorrelationFunction.h"


using namespace std;
using namespace Eigen;


VectorXd GetCorrelationVector(
  const Config &config, 
  BasisSet &atomicBasis, 
  const vector<size_t> &canonicalSortedLatticeIds, 
  const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters);

#endif // LMC_CE_INCLUDE_CORRELATIONVECTOR_H_

