#ifndef LMC_LMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_
#define LMC_LMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_
#include <random>
#include <omp.h>
#include "JumpEvent.h"
#include "TiledSupercell.h"
#include "KineticMcAbstract.h"
#include "VacancyMigrationPredictorTLMC.h"

namespace mc
{
  class KineticMcFirstOmp : public KineticMcFirstAbstract
  {
  public:
    KineticMcFirstOmp(TiledSupercell tiledSupercell,
                      const unsigned long long int logDumpSteps,
                      const unsigned long long int configDumpSteps,
                      const unsigned long long int maximumSteps,
                      const unsigned long long int thermodynamicAveragingSteps,
                      const unsigned long long int restartSteps,
                      const double restartEnergy,
                      const double restartTime,
                      const double temperature,
                      VacancyMigrationPredictorTLMC &vacancyMigrationPredictor,
                      const string &timeTemperatureFilename,
                      unique_ptr<RateCorrector> &rateCorrector,
                      const Eigen::RowVector3d &vacancyTrajectory);
    ~KineticMcFirstOmp() override;

  protected:
    void BuildEventList() override;
    double CalculateTime() override;
  };
} // mc

#endif // LMC_LMC_MC_INCLUDE_KINETICMCFIRSTOMP_H_
