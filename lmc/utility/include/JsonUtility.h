#ifndef _LMC_UTILITY_INCLUDE_JSONUTILITY_H_
#define _LMC_UTILITY_INCLUDE_JSONUTILITY_H_

#include <vector>
#include <omp.h>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
using namespace std;
using namespace Eigen;

/*! \brief Reads the parameters from trained data from json file
 *  \param json_filename Path to the JSON file.
 *  \param json_key Key under which values are stored.
 *  \return         Pair which containing the values for given json key.
 */
pair<VectorXd, double>
ReadParametersFromJson(const string &json_filename,
                       const string &json_key);

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

#endif // _LMC_UTILITY_JsonUtility_H_
