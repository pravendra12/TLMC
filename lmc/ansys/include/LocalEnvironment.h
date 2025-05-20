#ifndef LMC_ANSYS_INCLUDE_LOCALENVIRONMENT_H_
#define LMC_ANSYS_INCLUDE_LOCALENVIRONMENT_H_

#include "Config.h"
#include "ClusterExpansion.h"

using namespace std;
using namespace Eigen;

class LocalEnvironment
{
public:
  LocalEnvironment(
      Config config,
      size_t latticeId,
      const vector<double> &cutoffs);

  Config GetLocalConfig();

  VectorXd GetLocalConfigEncoding(
    const size_t maxClusterSize, 
    const size_t maxBondOrder);


private:
  const Config config_;
};

#endif // LMC_ANSYS_INCLUDE_LOCALENVIRONMENT_H_
