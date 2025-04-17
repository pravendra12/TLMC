#ifndef LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONBARRIERPREDICTOR_H_
#define LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONBARRIERPREDICTOR_H_

#include "Config.h"
#include "JsonUtility.h"
#include "Symmetry.h"
#include "EncodingUtility.h"
#include <omp.h>
#include <mutex>
#include <unordered_map>
#include <Eigen/Dense>
#include <boost/functional/hash.hpp>

using namespace std;
using namespace Eigen;
using PairMap = unordered_map<pair<size_t, size_t>,
                              vector<size_t>,
                              boost::hash<pair<size_t, size_t>>>;

class VacancyMigrationBarrierPredictor
{
public:
  // Constructor
  VacancyMigrationBarrierPredictor(const Config &config,
                                   const set<Element> &elementSet,
                                   const string &predictorFilename);

  /*! \brief Computes the Vacancy Migration Barrier.
   *  \param config A reference to the Config object containing lattice
                    configurations and relative positions.
   *  \param latticeIdJumpPair    Lattice Id jump pair.
  */
  [[nodiscard]] double GetBarrier(
      const Config &config,
      const pair<size_t, size_t> &latticeIdJumpPair) const;

private:
  const pair<VectorXd, double> barrier_fitted_parameters_{};

  /// Adjusted Effective cluster interaction.
  const VectorXd adjusted_beta_barrier_{};

  /// Adjusted Intercept
  const double adjusted_intercept_barrier_{};

  /// One Hot Encoding Hashmap
  const unordered_map<string, RowVectorXd> oneHotEncodingMap_{};

  /// Equivalent site encoding based on 3 Fold Rotation
  const vector<vector<size_t>> encoding3FoldRotation_{};

  // Currently considers neighbours upto 2nd bond order for the jump pair
  // Neighbors upto maxBondOrder_ are currently considered for making the
  // barrier prediction
  static constexpr size_t maxBondOrderForBarrierPrediction_ = 2;

protected:
  // For efficieny store all the symmetrically sorted vector once in a map
  const PairMap symmetricallySortedVectorMap_{};
};

/*! \brief Computes the symmetrically sorted vector map
 *
 *  This function utilizes multithreading and mutex to compute the symmetrically
 *  sorted vector map.
 *
 *  \param config A reference to the Configuration object
 *  \param maxBondOrder The maximum bond order; neighbors up to this bond order
 *  are used to generate the symmetrically sorted vector map
 *  \param symmetricallySortedVectorMap A map that associates lattice ID pairs
 *  with their corresponding symmetrically sorted vectors
 */
PairMap ComputeSymmetricallySortedVectorMap(const Config &config,
                                            const size_t maxBondOrder);

#endif // LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_
