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
#include "Symmetry.h"
#include "GetOrbits.h"
#include "JsonUtility.h"
#include "LocalEnvironmentEncoder.h"

using namespace std;
using namespace Eigen;

using PairMap = unordered_map<pair<size_t, size_t>,
                              vector<size_t>,
                              boost::hash<pair<size_t, size_t>>>;

class KRAPredictor
{
public:
    /*! @brief KRAPredictor Constructor
        @param predictorFilename Predictor filename which contains the fitting coefficient
        @param referenceConfig   Reference config
        @param elementSet        Element Set
     */
    KRAPredictor(const string &predictorFilename,
                 const Config &referenceConfig,
                 const set<Element> &elementSet);

    /*! @brief Function to compute Kinetically resolved activation barrier
        @param config            Configuration used for computing KRA
        @param latticeIdJumpPair Lattice id jump pair
        @return Kinetically resolved activation barrier
     */
    [[nodiscard]] double GetKRA(
        const Config &config,
        const pair<size_t, size_t> &latticeIdJumpPair) const;

private:
    /*! @brief Fitting coefficients
     */
    const VectorXd betaKRA_W_;
    const double interceptKRA_W_;

    const VectorXd betaKRA_Ta_;
    const double interceptKRA_Ta_;

    /*! @brief Element set
     */
    const set<Element> elementSet_;

    /*! @brief KRA from predictor file
     */
    const size_t maxBondOrder_;
    const size_t maxBondOrderOfCluster_;
    const size_t maxClusterSize_;

    /*! @brief Equivalent sites encoding under 3 Bar symmetry
     */
    const vector<vector<size_t>> equivalentSiteEncoding_;

    /*! @brief Canonical latticeIdPair Map to ssVector
     */
    const PairMap latticePairToSSVectorMap_;

    /*! @brief Atomic Basis
     */
    const string basisType_ = "Chebyshev";
};

/*!
 * @brief Returns a map of symmetrically sorted lattice ID pairs under 3̅ (3-bar) symmetry
 *        for all canonical lattice ID pairs.
 *
 * @param config        Configuration object used to define the lattice.
 * @param maxBondOrder  Maximum bond order used for neighborhood determination.
 *
 * @return PairMap containing symmetrically sorted lattice ID pairs.
 */
PairMap GetSymmetricallySortedLatticePairMap(
    const Config &config,
    const size_t maxBondOrder);

#endif // LMC_PRED_INCLUDE_KRAPREDICTOR_H_
