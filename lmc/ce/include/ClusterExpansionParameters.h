/*******************************************************************************
 * Copyright (c) 2025. All rights reserved.
 * @Author: Pravendra Patel
 * @Date:    2025-08-19
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-08-19
 *******************************************************************************/

#ifndef LMC_CE_INCLUDE_CLUSTEREXPANSIONPARAMETERS_H_
#define LMC_CE_INCLUDE_CLUSTEREXPANSIONPARAMETERS_H_

#include <vector>
#include <omp.h>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <iostream>
#include <nlohmann/json.hpp>
#include <set>
#include <map>
#include <iomanip>
#include "Element.hpp"
#include "JsonUtility.h"
#include <boost/functional/hash.hpp>

using json = nlohmann::json;
using namespace std;
using namespace Eigen;

/**
 * @class ClusterExpansionParameters
 * @brief Reads and stores cluster expansion and Kinetically Resolved Activation
 *        Barrier (KRA) parameters from a JSON coefficients file.
 *
 * This class encapsulates all the parameters required for a cluster expansion (CE) model
 * and KRA calculations, including:
 *   - Cluster expansion parameters (CE)
 *   - Kinetically resolved activation barrier parameters (KRA)
 *   - Kinetic effective coefficients of interactions (KECIs)
 *
 * It parses a JSON coefficients file with the following format:
 *
 * {
 *   "maxBondOrderKRA" : 2,
 *   "maxBondOrderOfClusterKRA" : 3,
 *   "maxClusterSizeKRA" : 3,
 *   "maxBondOrderCE" : 3,
 *   "maxClusterSizeCE" : 3,
 *   "basisType": "Occupation",
 *   "ce" : {
 *       "elementSet" : ["A", "B", "X"],
 *       "eci" : [...]
 *   },
 *   "kra" : {
 *       "elementSet" : ["A", "B"],
 *       "keci" : {
 *           "A" : [],
 *           "B" : []
 *       }
 *   }
 * }
 *
 * @note The class optionally provides a debug function that prints all parameters to the console.
 */
class ClusterExpansionParameters
{
public:
  /**
   * @brief Constructor that loads and parses the coefficients JSON file.
   *
   * @param coefficientFilename Path to the JSON coefficients file.
   * @param debug Optional flag to enable debug printing of all parameters after loading.
   */
  ClusterExpansionParameters(const string &coefficientFilename, const bool debug = false);

  string GetCoefficientFile() const;

  /**
   * @brief Returns the basis type used for the cluster expansion.
   * @return Basis type string (e.g., "Occupation" or "Chebyshev").
   */
  string GetBasisType() const;

  size_t GetMaxBondOrder(const string &jsonKey) const;

  size_t GetMaxClusterSize(const string &jsonKey) const;

  // can be use with kra / ce / lvfe
  set<Element> GetElementSet(const string &jsonKey) const;

  /**
   * @brief Returns the Effective Cluster Interaction (ECI) coefficients for the CE model.
   * @param ceType `ce` or `symCE`
   * @return VectorXd containing all ECIs.
   */
  vector<double> GetECIs(const string &ceType) const;

  // ------------------- Symmeteric Cluster Expansion -------------------
  /*

  symCE :
  {
    clusterCutoffs : [],
    allowedElements : [],
    chemicalPotentials :
    {
      A : muA,
      B : muB
    }
    ecis : []
  }

  */

  // {pair, triplet, quadraplet}
  vector<double> GetClusterCutoffs() const;

  // Allowed elements
  vector<string> GetAllowedElements() const;

  unordered_map<string, double> GetChemicalPotentialsMap() const;

  // GetECIs will remain only the jsonKey will be now taken as parameter

  // ------------------- Kinetically Resolved Activation Barrier (KRA) -------------------
  Vector3d GetReferenceJumpDirection() const;

  /**
   * @brief Returns the Kinetic Effective Coefficients of Interactions (KECIs) for each element.
   *
   * @return An unordered_map mapping each Element to its corresponding VectorXd of KECI values.
   */
  unordered_map<Element, VectorXd, boost::hash<Element>> GetKECIsMap() const;

  // unordered_map<Element, VectorXd, boost::hash<Element>> GetKECIs() const;
  VectorXd GetKECIs() const;

private:
  /**
   * @brief Debug function that prints all member function outputs to the console.
   *
   * Can be used to verify that all parameters were correctly loaded.
   */
  void DebugAllFunctions();

  json allParameters_; /**< JSON object storing all parameters */

  const string predictorFilename_;
};

#endif // LMC_CE_INCLUDE_CLUSTEREXPANSIONPARAMETERS_H_
