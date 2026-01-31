#pragma once

#include "LVFEPredictor.h"
#include "SymmetricCEPredictor.h"
#include "Config.h"

using namespace std;

// Computes Vacancy Formation Energy for central site in the config
map<Element, double> ComputeVFE(
    Config &config,
    const unordered_map<string, double> &chemicalPotentialMap,
    SymmetricCEPredictor &symCEPredictor,
    LVFEPredictor *lvfePredictor = nullptr);
