/**
 * @file CorrelationFunction.cpp
 * @brief Defines the Correlation Function for lattice Monte Carlo simulations.
 *
 * This file implements methods and members of the AtomBasis class, which
 * handles atomic basis configurations in simulations. It provides tools
 * for representing and manipulating atomic arrangements efficiently.
 *
 * @version 0.1
 * @date 2025-05-04
 * @copyright Copyright (c) 2025
 */

#include "CorrelationFunction.h"

RowVectorXd GetTensorProduct(const vector<RowVectorXd> &basisVector,
                             bool isSymmetric)
{
  if (basisVector.empty())
    return RowVectorXd();

  // If only one basis vector, return it directly
  if (basisVector.size() == 1)
    return basisVector[0];

  // For symmetric product: average over all permutations
  // Cluster which are formed by sites which are equivalent due to symmetry
  // W-Ta and Ta-W will be equivalent if they are formed by sites which are
  // equivalent.
  //
  // W-Ta : [0 0 1 0]
  // Ta-W : [0 1 0 0]
  //
  // Due to symmetry basis vector for W-Ta will be
  // 1/2 * (basis(W)*basis(Ta) + basis(Ta)*basis(W))
  // which will result into [0 0.5 0.5 0]
  if (isSymmetric)
  {
    vector<int> indices(basisVector.size());
    for (size_t i = 0; i < indices.size(); ++i)
      indices[i] = i;

    std::sort(indices.begin(), indices.end());

    RowVectorXd sum;
    int count = 0;

    do
    {
      vector<RowVectorXd> permuted;
      for (int idx : indices)
        permuted.push_back(basisVector[idx]);

      RowVectorXd prod = GetTensorProduct(permuted, false);

      if (sum.size() == 0)
        sum = RowVectorXd::Zero(prod.size());

      sum += prod;
      ++count;

    } while (std::next_permutation(indices.begin(), indices.end()));

    return sum / static_cast<double>(count);
  }

  // Standard (non-symmetric) tensor product
  // Cluster which are formed by sites which are not equivalent.
  // W-Ta and Ta-W will not be equivalent
  // W-Ta : [0 0 1 0]
  // Ta-W : [0 1 0 1]
  RowVectorXd result = basisVector[0];

  for (size_t i = 1; i < basisVector.size(); ++i)
  {
    const RowVectorXd &vec = basisVector[i];
    RowVectorXd temp(result.size() * vec.size());

    for (int j = 0; j < result.size(); ++j)
    {
      temp.segment(j * vec.size(), vec.size()) = result(j) * vec;
    }

    result = temp;
  }

  return result;
}

RowVectorXd GetTensorProduct(const RowVectorXd &basisVector1,
                             const RowVectorXd &basisVector2)
{
  // Outer product
  MatrixXd outerProduct = basisVector1.transpose() * basisVector2;

  // Flatten row-wise
  return Map<RowVectorXd>(outerProduct.data(), outerProduct.size());
}

// Return correlation function for an orbit
RowVectorXd GetCorrelationFunction(const Config &config,
                                   const set<Element> &elementSet,
                                   const string &basisType,
                                   const vector<vector<size_t>> &orbitVector,
                                   const bool &isClusterSymmetric)

{

  // Φ
  RowVectorXd corrFunction;
  bool isCorrFunctionResized = false;

  // Will be reused as many times as the function is recalled
  static AtomBasis atomicBasis(elementSet, basisType);

  // Number of clusters in the orbit
  int numClusters = 0;

  for (auto cluster : orbitVector)
  {
    vector<RowVectorXd> atomBasisVector;

    vector<string> elementCluster;

    // Retrieve the basis vector for the cluster
    // Iterate over the encoded cluster to extract the basis for each element
    // and compute the tensor product of the basis vectors
    for (auto latticeId : cluster)
    {
      auto element = config.GetElementOfLattice(latticeId);
      auto elementString = element.GetElementString();

      elementCluster.emplace_back(elementString);

      // Get the basis vector
      RowVectorXd atomBasis = atomicBasis.GetCachedAtomBasis(element);

      atomBasisVector.emplace_back(atomBasis);

      // cout << element.GetElementString() << latticeId << "( " << idx << " )" << " " << atomBasis << endl;
    }

    // Key for caching the tensor product
    string elementClusterString;

    if (isClusterSymmetric)
    {
      sort(elementCluster.begin(), elementCluster.end());
      for (const auto &element : elementCluster)
      {
        elementClusterString += element;
      }
    }
    else
    {
      // A-B and B-A will be different
      for (const auto &element : elementCluster)
      {
        elementClusterString += element + "-";
      }

      if (!elementClusterString.empty())
        elementClusterString.pop_back();
    }

    RowVectorXd clusterBasisVector = atomicBasis.GetCachedTensorProduct(
        elementClusterString,
        atomBasisVector,
        isClusterSymmetric);

    if (!isCorrFunctionResized)
    {
      corrFunction.resize(clusterBasisVector.size());
      corrFunction.setZero();
      isCorrFunctionResized = true;
    }

    // Add the cluster basis vector to the correlation function
    corrFunction += clusterBasisVector;
    numClusters++;
  }

  // Normalize the correlation function by the number of clusters in an orbit
  corrFunction /= numClusters;

  return corrFunction;
}

/*
RowVectorXd GetCorrelationFunction(const Config &config,
                                   const set<Element> &elementSet,
                                   const string &basisType,
                                   const vector<vector<size_t>> &orbitVector,
                                   const bool &isClusterSymmetric)

{

  // Φ
  RowVectorXd corrFunction;
  bool isCorrFunctionResized = false;

  // Number of clusters in the orbit
  int numClusters = 0;

  for (auto cluster : orbitVector)
  {
    vector<RowVectorXd> atomBasisVector;

    // cout << "-------------------" << endl;
    vector<string> elementCluster;

    // Retrieve the basis vector for the cluster
    // Iterate over the encoded cluster to extract the basis for each element
    // and compute the tensor product of the basis vectors
    for (auto latticeId : cluster)
    {
      // auto latticeId = symmetricSortedVector[idx];
      auto element = config.GetElementOfLattice(latticeId);

      elementCluster.emplace_back(element.GetElementString());

      RowVectorXd atomBasis = GetAtomBasis(element,
                                           elementSet,
                                           basisType);

      atomBasisVector.emplace_back(atomBasis);

      // cout << element.GetElementString() << latticeId << "( " << idx << " )" << " " << atomBasis << endl;
    }

    // cout << "Size of atomBasisVector: " << atomBasisVector.size() << endl;
    // cout << elementCluster << endl;

    RowVectorXd clusterBasisVector = GetTensorProduct(atomBasisVector,
                                                      isClusterSymmetric);

        // cout << clusterBasisVector << endl;

    if (!isCorrFunctionResized)
    {
      corrFunction.resize(clusterBasisVector.size());
      corrFunction.setZero();
      isCorrFunctionResized = true;
    }

    // Add the cluster basis vector to the correlation function
    corrFunction += clusterBasisVector;
    numClusters++;
  }

  // Normalize the correlation function by the number of clusters in an orbit
  corrFunction /= numClusters;
  // cout << "Normalized Correlation Function: " << corrFunction << endl;

  return corrFunction;
}
*/
