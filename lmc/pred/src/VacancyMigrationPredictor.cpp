#include "VacancyMigrationPredictor.h"

VacancyMigrationPredictor::VacancyMigrationPredictor(
    const string &predictorFilename,
    const Config &referenceConfig,
    const Config &supercellConfig,
    const set<Element> &elementSet,
    size_t maxClusterSize,
    size_t maxBondOrder) : eKRAPredictor(predictorFilename,
                                         referenceConfig,
                                         elementSet),
                           energyChangePredictor_(predictorFilename,
                                                  referenceConfig,
                                                  supercellConfig,
                                                  elementSet)
{
}

// (barrier, dE)
pair<double, double> VacancyMigrationPredictor::GetBarrierAndDeltaE(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair) const
{
  double eKRA = eKRAPredictor.GetKRA(config,
                                     latticeIdJumpPair);

  double dE = energyChangePredictor_.GetDeThreadSafe(config,
                                                     latticeIdJumpPair);

  // EKRA = Ea - 1/2*dE

  double barrier = eKRA + (dE / 2);

  return pair<double, double>(barrier, dE);
}

static std::unordered_map<ClusterType, size_t, boost::hash<ClusterType>> ConvertSetToHashMap(
    const std::set<ClusterType> &cluster_type_set)
{
  std::unordered_map<ClusterType, size_t, boost::hash<ClusterType>> cluster_type_count;
  for (const auto &cluster_type : cluster_type_set)
  {
    cluster_type_count[cluster_type] = 0;
  }
  return cluster_type_count;
}
