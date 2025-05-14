#ifndef _LMC_UTILITY_INCLUDE_JSONUTILITY_H_
#define _LMC_UTILITY_INCLUDE_JSONUTILITY_H_

#include <vector>
#include <omp.h>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <iostream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
using namespace std;
using namespace Eigen;

/*! @brief Reads the parameter from the json file.

    Even if there is only one value corresponding to the subkey in the parameter
    a vector will be returned.

    @param jsonFilename
    @param jsonKey
    @param jsonSubKey
    @return Vector containing the values.
*/
VectorXd ReadParametersFromJson(const string &jsonFilename,
                                const string &jsonKey,
                                const string &jsonSubKey);

size_t ReadParameterFromJson(const string &jsonFilename,
                             const string &jsonKey);

#endif // _LMC_UTILITY_JsonUtility_H_
