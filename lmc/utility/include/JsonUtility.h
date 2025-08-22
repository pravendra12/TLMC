#ifndef _LMC_UTILITY_INCLUDE_JSONUTILITY_H_
#define _LMC_UTILITY_INCLUDE_JSONUTILITY_H_

#include <vector>
#include <omp.h>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <iostream>
#include <nlohmann/json.hpp>
#include "Element.hpp"
#include <set>
#include <map>

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

vector<string> ReadStringParametersFromJson(const string &jsonFilename,
                                            const string &jsonKey,
                                            const string &jsonSubKey);

size_t ReadParameterFromJson(const string &jsonFilename,
                             const string &jsonKey);

string ReadBasisType(const string &jsonFilename,
                     const string &jsonKey);

void ReadKRAParametersFromJson(
    const string &jsonFilename);

template <typename T>
T ReadParameterFromJson(
    const json &parameters,
    const string &jsonKey);

template <typename T>
vector<T> ReadParametersFromJson(
    const json &allParameters,
    const string &jsonKey,
    const string &jsonSubKey);

#include "JsonUtility.inl"

#endif // _LMC_UTILITY_JsonUtility_H_
