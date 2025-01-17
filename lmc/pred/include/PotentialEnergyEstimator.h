/**************************************************************************************************
 * Copyright (c) 2022-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 11/15/22 12:36 PM                                                                       *
 * @Last Modified by: pravendra12                                                                 *
 * @Last Modified time: 11/30/24 5:00 PM                                                          *
 **************************************************************************************************/

/*! \file  PotentialEnergyEstimator.h
 *  \brief File for the PotentialEnergyEstimator class definition.
 */
#ifndef LMC_PRED_INCLUDE_POTENTIALENERGYESTIMATOR_H_
#define LMC_PRED_INCLUDE_POTENTIALENERGYESTIMATOR_H_

#include <set>
#include <numeric>
#include <algorithm>

#include "eigen3/Eigen/Dense"
#include "ClusterExpansion.h"

/*! \brief Class for defining cluster expansion Hamiltonian.
 */
class PotentialEnergyEstimator {
 public:
  PotentialEnergyEstimator(const std::string &predictor_filename,
                           const Config &reference_config,
                           const Config &supercell_config,
                           const std::set<Element> &element_set,
                           size_t max_cluster_size,
                           size_t max_bond_order);
  ~PotentialEnergyEstimator();

  /*! \brief Get the encode vector of the configuration, which is the number of 
   *         appearance of each cluster types plus void cluster.
   *  \param config   The configuration the code works on
   *  \return         The encode vector
   */  
  [[nodiscard]] Eigen::VectorXd GetEncodeVector(const Config &config) const;
  
  /*! \brief Get the encoded vector of the configuration, representing the count
   *         of each cluster type, including void clusters.
   *  \param config           The configuration to analyze.
   *  \param lattice_cluster  Vector containing lattice site IDs for the cluster.
   *  \return                 A vector representing the encoded cluster type counts.
   */
  [[nodiscard]] Eigen::VectorXd GetEncodeVectorOfCluster(const Config &config, 
                                                        std::vector<size_t> cluster) const;

  /*! \brief Get the Energy of the configuration
   *  \param config           The configuration to analyze.
   *  \return                 Energy of the configuration.
   */                                                      
  [[nodiscard]] double GetEnergy(const Config &config) const;

  /*! \brief Get the energy of the cluster.
   *  \param config           The configuration to analyze.
   *  \param cluster          Vector containing lattice site IDs for the cluster.
   *  \return                 Energy of the cluster.
   */
  [[nodiscard]] double GetEnergyOfCluster(const Config &config,
                                          const std::vector<size_t> &cluster) const;
  
  /*! \brief Get the energy change due to atom jump.
   *  \param config           The configuration to analyze.
   *  \param lattice_id_pair  Lattice Id jump pair.
   *  \return                 Change in energy due to atom jump.
   */
  [[nodiscard]] double GetDe(Config &config, const std::pair<size_t, size_t> &lattice_id_pair) const;
  [[nodiscard]] std::map<Element, double> GetChemicalPotential(Element solvent_element) const;

  /// Effective cluster interaction.
  const Eigen::VectorXd effective_cluster_interaction_{};

  
 private:
 
  /// Element set
  const std::set<Element> element_set_{};

  /// Set that contains all available cluster types
  const std::set<ClusterType> initialized_cluster_type_set_{};   

  /// Maps each lattice cluster type to its count. 
  /// Used mainly for configuration used for training Cluster Expansion Model.
  const std::unordered_map<LatticeClusterType, size_t, 
                boost::hash<LatticeClusterType>> lattice_cluster_type_count_{};
  
  /// Maximum Cluster Size
  const size_t max_cluster_size_{};

  /// Maximum Bond Order
  const size_t max_bond_order_{};
};

#endif //LMC_PRED_INCLUDE_POTENTIALENERGYESTIMATOR_H_