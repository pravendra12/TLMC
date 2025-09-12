/*******************************************************************************
 * Copyright (c) 2025. All rights reserved.
 * @Author: Pravendra Patel
 * @Date:    2025-06-01
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! @file KRAPredictor.h
    @brief File contains declaration of Kinetically Resolved Activation Barrier
*/

#ifndef LMC_PRED_INCLUDE_KRAPREDICTOR_H_
#define LMC_PRED_INCLUDE_KRAPREDICTOR_H_

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

using namespace std;
using namespace Eigen;

class KRAPredictor
{
public:
  /*! @brief KRAPredictor Constructor
      @param predictorFilename Predictor filename which contains the fitting coefficient
      @param referenceConfig   Reference config
      @param elementSet        Element Set
   */

  KRAPredictor(
      const ClusterExpansionParameters &ceParams,
      const Config &config);

  /*! @brief Function to compute Kinetically resolved activation barrier
      @param config            Configuration used for computing KRA
      @param latticeIdJumpPair Lattice id jump pair
      @return Kinetically resolved activation barrier
   */
  [[nodiscard]] double GetKRA(
      const Config &config,
      const pair<size_t, size_t> &latticeIdJumpPair) const;

private:
 
  // Assign a pair a index which can be used to access the symmetrically sorted sites
  // around that lattice id pair
  using PairIndexMap = unordered_map<pair<size_t, size_t>, size_t, boost::hash<pair<size_t, size_t>>>;

  const size_t maxBondOrder_{};
  const size_t maxClusterSize_{};

  // Keep this as mutable allowing const function to make changes to it
  mutable BasisSet atomicBasis_;

  const Vector3d referenceJumpDirection_{};

  const unordered_map<size_t, RowVector3d> canonicalReferenceMap_{};

  // Contains space group symmetry operation
  const vector<pair<Matrix3d, Vector3d>> symmetryOperations_{};

  /*! @brief Equivalent Clusters Encoding
   */
  const vector<pair<vector<vector<size_t>>,
                    LatticeClusterType>>
      encodedOrbitsForPair_{};

  const unordered_map<Element, VectorXd, boost::hash<Element>> KECIsMap_{};

  vector<vector<size_t>> symmetricallySortedLatticeIdsVectorMap_{};
  PairIndexMap pairToIndexHashMap_{};


  void GetSymmetricallySortedLatticeIdsVectorMap(
    const Config &config);

  void PrintKRAPredictorInfo() const;

};

#endif // LMC_PRED_INCLUDE_KRAPREDICTOR_H_
