#ifndef LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORLRU_H_
#define LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORLRU_H_

#include "VacancyMigrationPredictor.h"
#include "LruCache.hpp"

using namespace std;

class VacancyMigrationPredictorLru
{
public:
VacancyMigrationPredictorLru(const Config &config,
                                      const set<Element> &elementSet,
                                      const string &predictorFilename,
                                      size_t cacheSize);

  // ~VacancyMigrationPredictorLru() override;

  // [[nodiscard]] pair<double, double>
  // GetBarrierAndDiffFromLatticeIdPair(const Config &config,
  //                                    const pair<size_t, size_t> &latticeIdJumpPair) const;

private:
  [[nodiscard]] size_t GetHashFromConfigAndLatticeIdPair(const Config &config,
                                                         const pair<size_t, size_t> &latticeIdJumpPair) const;

  mutable LruCache<size_t, std::pair<double, double>> lruCache_;
};

#endif //LMC_PRED_INCLUDE_VACANCYMIGRATIONPREDICTORLRU_H_
