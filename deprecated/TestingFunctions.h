#ifndef LMC_UTILITY_INCLUDE_TESTINGFUNCTIONS_H_
#define LMC_UTILITY_INCLUDE_TESTINGFUNCTIONS_H_

#include "Config.h"
#include "EnergyPredictor.h"
#include "PrintUtility.h"


using namespace std;

void TestEnergyChangePredictor(
  Config &config, 
  const string &predictorFilename);

#endif // LMC_UTILITY_INCLUDE_TESTINGFUNCTIONS_H_

