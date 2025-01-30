#ifndef LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_
#define LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_

// #include "LruCache.hpp"
#include "Config.h"
#include "PotentialEnergyEstimator.h"



class VacancyMigrationPredictor {
 public:
  // Constructor
  VacancyMigrationPredictor(const std::string &predictor_filename);

//  [[nodiscard]] std::pair<double, double>
//  GetBarrierAndDiffFromLatticeIdPair(Config &config,
//                                     const std::pair<size_t, size_t> &lattice_id_jump_pair) const;

  [[nodiscard]] std::array<double, 3>
  GetBarrierAndDiffFromLatticeIdPair(
                 Config &config,
                 const std::pair<size_t, size_t> &lattice_id_jump_pair) const;
  
  /*! \brief Computes the Vacancy Migration Barrier.
   *  \param config                  A reference to the Config object containing 
   *                                 lattice configurations and relative positions.
   *  \param lattice_id_jump_pair    Lattice Id jump pair.
  */
  [[nodiscard]] double GetBarrier(const Config &config, 
                                  const std::pair<size_t, size_t> &lattice_id_jump_pair) const;

  double GetBarrierNew(const Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair) const;

  double GetDiffNew(Config &config, const std::pair<size_t, size_t> lattice_id_jump_pair) const;

  /*! \brief Computes the driving force using the Df of the initial and final vacancy structure.
   *  \param config                  A reference to the Config object containing 
   *                                 lattice configurations and relative positions.
   *  \param lattice_id_jump_pair    Lattice Id jump pair.
  */
  [[nodiscard]] double GetDiff(Config &config, 
                      const std::pair<size_t, size_t> lattice_id_jump_pair) const;

  

  private:

  // For now these values are being written randomly.
  // But later these will read from predictor file.

  /// Slope for predicting driving force.
  double slope_E_diff_{}; 
  /// Intercept for predicting driving force.
  double intercept_E_diff_{}; 
  /// Slope for DbE.
  double slope_DbE_{};
  /// Slope for DbA.
  double slope_DbA_ ={};
  /// Intercept for barrier prediction.
  double intercept_barrier_{};

//   mutable LruCache<size_t, std::pair<double, double>> lru_cache_;

}; 


/*! \brief Computes the contribution due to electronic interaction of an atom.
 *  \param element    Element for which xSv to be computed.
 *  \return           The xSv of the migrating atom.
*/
double GetxSv(const Element &element);

/*!
 * \brief Computes the Df for the neighboring lattice sites.
 * \param config      A reference to the Config object containing lattice 
 *                    configurations and relative positions.
 * \param lattice_id  The lattice ID for which Df to be calculated.
 * \return            The Df value for the given lattice ID.
 */
double GetDf(const Config &config, const size_t lattice_id);

/*!
 * \brief Computes the geometric mean of a vector of positive values.
 * 
 * \param xSv_vector   A vector containing the values for which the geometric 
 *                     mean is to be calculated. All elements must be positive.
 * \param z            An optional scaling factor applied to the exponent. 
 *                     Default is 1.0.
 * \return             The geometric mean of the input values scaled by `z`.
 */
double GetGeometricMean(std::vector<double> &xSv_vector, double z);



#endif    //LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_
