#pragma once

#include "TiledSupercell.h"
#include "Config.h"
#include "Parameter.h"
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

// Path to TLMC folder where the atomIndicies.bin.gz files are present
// It will search for cmc_param.txt or kmc_param.txt to read the required 
// parameters
void ConvertAtomVectorsToConfigs(
  const string &pathToTLMCRun);