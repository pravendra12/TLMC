#ifndef LMC_CONSTANT_INCLUDE_SAVECLUSTERCONFIG_H_
#define LMC_CONSTANT_INCLUDE_SAVECLUSTERCONFIG_H_

#include "Config.h"
#include "ClusterExpansion.h"

using namespace std;

void SaveClustersConfig(
    const string outputPath,
    const Config &config,
    const vector<set<vector<size_t>>> &equivalentClusters);

#endif // LMC_CONSTANT_INCLUDE_SAVECLUSTERCONFIG_H_
