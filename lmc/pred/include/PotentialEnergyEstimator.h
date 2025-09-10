/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 6/14/22 12:36 PM
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! @file  PotentialEnergyEstimator.h
    @brief File contains declaration of Potential Energy Estimator Class
 */

#ifndef LMC_PRED_INCLUDE_POTENTIALENERGYESTIMATOR_H_
#define LMC_PRED_INCLUDE_POTENTIALENERGYESTIMATOR_H_

#include <set>
#include <numeric>
#include <algorithm>
#include <Eigen/Dense>
#include "ClusterExpansion.h"
#include "ClusterExpansionParameters.h"

using namespace std;
using namespace Eigen;

class PotentialEnergyEstimator
{
public:
    /*! @brief Constructor for Potential Energy Estimator
        @param ceParams          Contains Cluster Expansion Parameters
        @param referenceConfig   Configuration used for defining the object
        @param supercellConfig   Training supercell size config used for training
    */
    PotentialEnergyEstimator(
        const ClusterExpansionParameters &ceParams,
        const Config &referenceConfig,
        const Config &trainingConfig);

    /*! @brief Destructor for Potential Estimator
     */
    ~PotentialEnergyEstimator();

    /*! @brief Get the encode vector of the configuration, which is the number of
     *         appearance of each cluster types plus void cluster.
     *  @param config   The configuration the code works on
     *  @return         The encode vector for the entire configuration
     */
    [[nodiscard]] VectorXd GetEncodeVector(const Config &config) const;

    /*! @brief Get the encoded vector of the configuration, representing the count
     *         of each cluster type, including void clusters.
     *  @param config           The configuration to analyze.
     *  @param lattice_cluster  Vector containing lattice site IDs for the cluster.
     *  @return                 A vector representing the encoded cluster type counts.
     */
    [[nodiscard]] VectorXd GetEncodeVectorOfCluster(
        const Config &config,
        const std::vector<size_t> &cluster) const;

    /*! @brief Get the Energy of the configuration
     *  @param config           The configuration to analyze.
     *  @return                 Energy of the configuration.
     */
    [[nodiscard]] double GetEnergy(const Config &config) const;

    /*! @brief Get the energy of the cluster.
     *  @param config           The configuration to analyze.
     *  @param cluster          Vector containing lattice site IDs for the cluster.
     *  @return                 Energy of the cluster.
     */
    [[nodiscard]] double GetEnergyOfCluster(
        const Config &config,
        const vector<size_t> &cluster) const;

    // Take allowed sites lattice Id vectors
    // Need to change the name of function
    Eigen::VectorXd GetEncodeVectorWithinAllowedSites(const Config &config, const vector<size_t> &allowedSites) const;

    // Returns the energy of the cluster formed but allowed sites
    // Need to think of better name but GetEnergyOfCluster already there
    // allowedSites contains the lattice Ids
    double GetEnergyOfClusterWithinAllowedSites(const Config &config, const vector<size_t> &allowedSites) const;

    /*! @brief Get the energy change due to atom swap.
     *  @param config           The configuration to analyze.
     *  @param latticeIdPair    Lattice Id swap pair.
     *  @return                 Change in energy due to atom swap.
     */
    [[nodiscard]] double GetDeSwap(
        Config &config,
        const pair<size_t, size_t> &latticeIdPair) const;

    /*! @brief Get the energy change due to atom migration.
     *  @param config               The configuration to analyze.
     *  @param latticeIdJumpPair    Lattice Id jump pair.
     *  @return                     Change in energy due to atom migration.
     */
    [[nodiscard]] double GetDeMigration(
        const Config &config,
        const std::pair<size_t, size_t> &latticeIdJumpPair) const;

private:
    /*! @brief Maximum Cluster Size
     */
    const size_t maxClusterSize_{};

    /*! @brief Maximum Bond Order
     */
    const size_t maxBondOrder_{};

    /*! @brief Effective cluster interactions
     */
    const VectorXd betaCE_{};

    /*! @brief Element set
     */
    const set<Element> elementSet_{};

    /*! @brief Set that contains all available cluster types
     */
    const set<ClusterType> initializedClusterTypeSet_{};

    /*! @brief Map that contains all available cluster types
     */
    const unordered_map<ClusterType, size_t, boost::hash<ClusterType>> clusterTypeCountHashMap_;

    /*! @brief Maps each lattice cluster type to its count. From config of the
     *         same size which was used for training Cluster Expansion Model.
     */
    const unordered_map<LatticeClusterType, size_t,
                        boost::hash<LatticeClusterType>>
        latticeClusterTypeCount_{};
};

#endif // LMC_PRED_INCLUDE_POTENTIALENERGYESTIMATOR_H_
