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
#include "Config.h"
#include "SymmetrySpglib.h"
#include "ClusterExpansionParameters.h"
#include "CorrelationVector.h"

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
        const pair<size_t, size_t> &latticeIdJumpPair);

private:
    const size_t maxBondOrder_;
    const size_t maxBondOrderOfCluster_;
    const size_t maxClusterSize_;

    BasisSet atomicBasis_;

    const VectorXd kecis_;

    const Vector3d referenceJumpDirection_{1, 1, 1};

    const unordered_map<size_t, RowVector3d> canonicalReferenceMap_;

    // Contains space group symmetry operation
    const vector<pair<Matrix3d, Vector3d>> symmetryOperations_;

    /*! @brief Equivalent Clusters Encoding
     */
    const vector<pair<vector<vector<size_t>>,
                      LatticeClusterType>>
        equivalentClustersEncoding_;
};

#endif // LMC_PRED_INCLUDE_KRAPREDICTOR_H_
