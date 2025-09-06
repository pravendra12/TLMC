#ifndef LMC_SYMCE_INCLUDE_INTERFACEUTILITY_H_
#define LMC_SYMCE_INCLUDE_INTERFACEUTILITY_H_

#include "Structure.hpp"
#include "Config.h"
#include <pybind11/numpy.h>
#include <pybind11/embed.h> 

using namespace std;
using namespace Eigen;
namespace py = pybind11;

// Converts Structure to Config
Config ConvertStructureToConfig(const Structure &structure);

// Converts Config to Structure
Structure ConvertConfigToStructure(const Config &config); 



#endif // LMC_SYMCE_INCLUDE_INTERFACEUTILITY_H_

