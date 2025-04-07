#ifndef LMC_PRED_INCLUDE_POTENTIALENERGYESTIMATOR_H_
#define LMC_PRED_INCLUDE_POTENTIALENERGYESTIMATOR_H_

#include <set>
#include <numeric>
#include <algorithm>
#include <Eigen/Dense>
#include "ClusterExpansion.h"
#include "JsonUtility.h"
using namespace std;
using namespace Eigen;

/*!
 * \file PotentialEnergyEstimator.h
 * \brief File for the PotentialEnergyEstimator class definition.
 * 
 * \author zhucongx
 * \date 6th March 2025
 * \lastedited 6th March 2025 by pravendra12
 */

/*! \brief Class for defining Cluster Expansion Hamiltonian.
 */
class PotentialEnergyEstimator {
 public:
  PotentialEnergyEstimator(const string &predictor_filename,
                           const Config &reference_config,
                           const Config &supercell_config,
                           const set<Element> &element_set,
                           size_t max_cluster_size,
                           size_t max_bond_order);
  ~PotentialEnergyEstimator();

  /*! \brief Get the encode vector of the configuration, which is the number of 
   *         appearance of each cluster types plus void cluster.
   *  \param config   The configuration the code works on
   *  \return         The encode vector
   */  
  [[nodiscard]] VectorXd GetEncodeVector(const Config &config) const;
  
  /*! \brief Get the encoded vector of the configuration, representing the count
   *         of each cluster type, including void clusters.
   *  \param config           The configuration to analyze.
   *  \param lattice_cluster  Vector containing lattice site IDs for the cluster.
   *  \return                 A vector representing the encoded cluster type counts.
   */
  [[nodiscard]] VectorXd GetEncodeVectorOfCluster(const Config &config, 
                                                  const std::vector<size_t> &cluster) const;

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
                                          const vector<size_t> &cluster) const;
  
  /*! \brief Get the energy change due to atom jump.
   *  \param config           The configuration to analyze.
   *  \param lattice_id_pair  Lattice Id jump pair.
   *  \return                 Change in energy due to atom jump.
   */
  [[nodiscard]] double GetDe(Config &config, const pair<size_t, size_t> &lattice_id_pair) const;

  double GetDeThreadSafe(const Config &config, const std::pair<size_t, size_t> &lattice_id_pair) const;

  
  [[nodiscard]] map<Element, double> GetChemicalPotential(Element solvent_element) const;

  private:
  const pair<VectorXd, double> ce_fitted_parameters_{};
  /// Adjusted Effective cluster interaction.
  const VectorXd adjusted_beta_ce_{};

  /// @brief Adjusted Intercept
  const double adjusted_intercept_ce_;
 
  /// Element set
  const set<Element> element_set_{};

  /// Set that contains all available cluster types
  const set<ClusterType> initialized_cluster_type_set_{};   

  /// Maps each lattice cluster type to its count. 
  /// Used mainly for configuration used for training Cluster Expansion Model.
  const unordered_map<LatticeClusterType, size_t, 
                boost::hash<LatticeClusterType>> lattice_cluster_type_count_{};
  
  /// Maximum Cluster Size
  const size_t max_cluster_size_{};

  /// Maximum Bond Order
  const size_t max_bond_order_{};

  
};

  
#endif //LMC_PRED_INCLUDE_POTENTIALENERGYESTIMATOR_H_