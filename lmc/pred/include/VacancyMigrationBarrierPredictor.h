#ifndef LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONBARRIERPREDICTOR_H_
#define LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONBARRIERPREDICTOR_H_

#include "Config.h"
#include "JsonUtility.h"
#include "Symmetry.h"
#include "EncodingUtility.h"
#include <unordered_map>
#include <eigen3/Eigen/Dense>
using namespace std;

class VacancyMigrationBarrierPredictor {
 public:
 
  // Constructor
  VacancyMigrationBarrierPredictor(const Config &config, 
                                   const set<Element> &elementSet,
                                   const size_t &maxBondOrder,
                                   const string &predictorFilename);
  

  /*! \brief Computes the Vacancy Migration Barrier.
   *  \param config                  A reference to the Config object containing 
   *                                 lattice configurations and relative positions.
   *  \param latticeIdJumpPair    Lattice Id jump pair.
  */
  [[nodiscard]] double GetBarrier(
              const Config &config, 
              const pair<size_t, size_t> &latticeIdJumpPair) const;

  private:
    const pair<Eigen::VectorXd, double> barrier_fitted_parameters_{};
    
    /// Adjusted Effective cluster interaction.
    const VectorXd adjusted_beta_barrier_{};

    /// Adjusted Intercept
    const double adjusted_intercept_barrier_{};
    
    /// One Hot Encoding Hashmap
    const unordered_map<string, RowVectorXd> oneHotEncodingMap_{};

    /// Equivalent site encoding based on 3 Fold Rotation 
    const vector<vector<size_t>> encoding3FoldRotation_{};


    // Future
    // For efficieny
    // Store all the symmetrically sorted vector once in a map

    const size_t maxBondOrder_{};







    

  




  
}; 

#endif    //LMC_LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTOR_H_
