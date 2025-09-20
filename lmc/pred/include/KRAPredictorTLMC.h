/*******************************************************************************
 * Copyright (c) 2025. All rights reserved.
 * @Author: Pravendra Patel
 * @Date:    2025-06-01
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! @file KRAPredictorTLMC.h
    @brief File contains declaration of Kinetically Resolved Activation Barrier
*/

#ifndef LMC_PRED_INCLUDE_KRAPREDICTORTLMC_H_
#define LMC_PRED_INCLUDE_KRAPREDICTORTLMC_H_

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <Eigen/Dense>
#include <boost/functional/hash.hpp>
#include "Config.h"
#include "SymmetrySpglib.h"
#include "CorrelationVector.h"
#include "GetEquivalentClusters.h"
#include "ClusterExpansionParameters.h"
#include "TiledSupercell.h"
#include "LatticeSiteMapping.hpp"
#include "LatticeSiteEncodedMapping.hpp"

using namespace std;
using namespace Eigen;

class KRAPredictorTLMC
{
public:
  KRAPredictorTLMC(
      const ClusterExpansionParameters &ceParams,
      const TiledSupercell &tiledSupercell);

  /*! @brief Function to compute Kinetically resolved activation barrier
        @param config            Configuration used for computing KRA
        @param latticeIdJumpPair Lattice id jump pair
        @return Kinetically resolved activation barrier
     */
  [[nodiscard]] double GetKRA(
      const TiledSupercell &tiledSupercell,
      const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair) const;

private:
  // Assign a pair a index which can be used to access the symmetrically sorted sites
  // around that lattice id pair
  using PairIndexMap = unordered_map<pair<size_t, size_t>, size_t, boost::hash<pair<size_t, size_t>>>;

  const size_t maxBondOrder_{};
  const size_t maxClusterSize_{};

  const Vector3d referenceJumpDirection_{};
  
  // Keep this as mutable allowing const function to make changes to it
  mutable BasisSet atomicBasis_;


  const unordered_map<size_t, RowVector3d> canonicalReferenceMap_{};

  // Contains space group symmetry operation
  const vector<pair<Matrix3d, Vector3d>> symmetryOperations_{};

  /*! @brief Equivalent Clusters Encoding
   */
  const vector<pair<vector<vector<size_t>>, LatticeClusterType>>
      encodedOrbitsForPair_{};

  const unordered_map<Element, VectorXd, boost::hash<Element>> KECIsMap_{};

  vector<vector<NeighbourOfPair>> symmetricLatticeSiteEncodedMappings_{};
  PairIndexMap pairToIndexHashMap_{};

  void GetSymmetricallySortedLatticeIdsVectorMap(
      const TiledSupercell &tiledSupercell);

  void PrintKRAPredictorInfo() const;
};

#endif // LMC_PRED_INCLUDE_KRAPREDICTORTLMC_H_
