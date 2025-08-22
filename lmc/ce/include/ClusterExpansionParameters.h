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

  /**
   * @brief Returns the basis type used for the cluster expansion.
   * @return Basis type string (e.g., "Occupation").
   */
  string GetBasisType() const;

  // ------------------- Cluster Expansion -------------------

  /**
   * @brief Returns the maximum bond order considered in the cluster expansion (CE).
   * @return Maximum bond order for CE.
   */
  size_t GetMaxBondOrderCE() const;

  /**
   * @brief Returns the maximum cluster size used in the cluster expansion (CE).
   * @return Maximum cluster size for CE.
   */
  size_t GetMaxClusterSizeCE() const;

  /**
   * @brief Returns the set of elements included in the cluster expansion (CE) model.
   * @return Set of elements in CE.
   */
  set<Element> GetElementSetCE() const;

  /**
   * @brief Returns the Effective Cluster Interaction (ECI) coefficients for the CE model.
   * @return VectorXd containing all ECIs.
   */
  VectorXd GetECIs() const;

  // ------------------- Kinetically Resolved Activation Barrier (KRA) -------------------

  /**
   * @brief Returns the maximum bond order of clusters considered in the KRA model.
   * @return Maximum bond order of clusters for KRA.
   */
  size_t GetMaxBondOrderOfClusterKRA() const;

  /**
   * @brief Returns the maximum bond order considered for KRA calculations.
   * @return Maximum bond order for KRA.
   */
  size_t GetMaxBondOrderKRA() const;

  /**
   * @brief Returns the maximum cluster size considered in the KRA model.
   * @return Maximum cluster size for KRA.
   */
  size_t GetMaxClusterSizeKRA() const;

  /**
   * @brief Returns the set of elements included in the KRA model.
   * @return Set of elements in KRA.
   */
  set<Element> GetElementSetKRA() const;

  /**
   * @brief Returns the Kinetic Effective Coefficients of Interactions (KECIs) for each element.
   *
   * @return An unordered_map mapping each Element to its corresponding VectorXd of KECI values.
   */
  unordered_map<Element, VectorXd, boost::hash<Element>> GetKECIs() const;

private:
  /**
   * @brief Debug function that prints all member function outputs to the console.
   *
   * Can be used to verify that all parameters were correctly loaded.
   */
  void DebugAllFunctions();
  
  json allParameters_; /**< JSON object storing all parameters */
};

#endif // LMC_CE_INCLUDE_CLUSTEREXPANSIONPARAMETERS_H_
