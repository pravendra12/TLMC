
/**
 * @file CorrelationFunction.h
 * @brief Provides declarations for correlation function use in Cluster Expansion.
 *
 * This header declares functions for computing tensor products of row vectors,
 * including symmetric and standard tensor products. These functions are used
 * to represent and manipulate atomic arrangements in simulations.
 *
 * @author Pravendra
 * @version 0.1
 * @date 2025-05-04
 * @copyright Copyright (c) 2025
 */

#ifndef LMC_CE_INCLUDE_CORRELATIONFUNCTION_H_
#define LMC_CE_INCLUDE_CORRELATIONFUNCTION_H_

#include <vector>
#include <Eigen/Dense>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <utility>
#include "AtomBasis.h"
#include "Config.h"
#include "LatticeCluster.hpp"

using namespace std;
using namespace Eigen;

/**
 * @brief Computes the tensor product of a set of row vectors.
 *
 * This function calculates the tensor product of the given basis vectors.
 * If `isSymmetric` is true, it computes the symmetric tensor product by
 * averaging over all permutations of the input vectors. Otherwise, it
 * computes the standard tensor product. The function efficiently handles
 * edge cases such as empty input or a single vector.
 *
 * Symmetric Tensor Product:
 * - Clusters are formed by sites equivalent due to symmetry.
 * - Example:
 *   - W-Ta and Ta-W are equivalent if formed by equivalent sites.
 *   - W-Ta : [0 0 1 0]
 *   - Ta-W : [0 1 0 0]
 *   - Symmetric basis vector for W-Ta:
 *     - 1/2 * (basis(W)*basis(Ta) + basis(Ta)*basis(W))
 *     - Result: [0 0.5 0.5 0]
 *
 * Standard (Non-Symmetric) Tensor Product:
 * - Clusters are formed by sites that are not equivalent.
 * - Example:
 *   - W-Ta and Ta-W are not equivalent.
 *   - W-Ta : [0 0 1 0]
 *   - Ta-W : [0 1 0 1]
 *
 * @param basisVector A vector of RowVectorXd representing the basis vectors.
 * @param isSymmetric A boolean flag to compute the symmetric tensor product.
 * @return RowVectorXd The resulting tensor product as a row vector.
 */
RowVectorXd GetTensorProduct(const vector<RowVectorXd> &basisVector,
                             bool isSymmetric);

/**
 * @brief Computes the tensor product of two row vectors.
 *
 * This function calculates the outer product of two row vectors and
 * flattens the result into a single row vector.
 *
 * @param basisVector1 The first row vector.
 * @param basisVector2 The second row vector.
 * @return RowVectorXd The resulting tensor product as a row vector.
 */
RowVectorXd GetTensorProduct(const RowVectorXd &basisVector1,
                             const RowVectorXd &basisVector2);

/**
 * @brief Computes the correlation function for a given orbit of clusters in a configuration.
 *
 * This function calculates the correlation function by taking the tensor product
 * of the basis vectors corresponding to the clusters in the specified orbit. It
 * supports both symmetric and non-symmetric tensor products based on the provided
 * basis type. The result is a row vector that represents the correlation function
 * for the given orbit in the configuration.
 *
 * @param config The configuration object representing the atomic arrangement.
 * @param elementSet A set of elements involved in the configuration.
 * @param basisType A string specifying the type of basis to use (e.g., symmetric or non-symmetric).
 * @param orbitVector A vector of vectors encoding the orbit of clusters.
 * @param isClusterSymmetric A boolean flag indicating whether the cluster is symmetric.
 * @param symmetricSortedVector A vector of sorted lattice Id's.
 * @return RowVectorXd A row vector representing the computed correlation function.
 */
RowVectorXd GetCorrelationFunction(const Config &config,
                                   const set<Element> &elementSet,
                                   const string &basisType,
                                   const vector<vector<size_t>> &orbitVector,
                                   const bool &isClusterSymmetric);

#endif // LMC_CE_INCLUDE_CORRELATIONFUNCTION_H_