#include "JsonUtility.h"

VectorXd ReadParametersFromJson(const string &jsonFilename,
                                const string &jsonKey,
                                const string &jsonSubKey)
{
  ifstream ifs(jsonFilename, ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  VectorXd parameterVector;

  for (const auto &[key, parameters] : all_parameters.items())
  {
    if (key == jsonKey)
    {
      for (const auto &[subKey, subParams] : parameters.items())
      {
        if (subKey == jsonSubKey)
        {
          vector<double> parameterArray;

          for (const auto &param : subParams)
          {
            parameterArray.push_back(param.get<double>());
          }

          // Convert the vector to VectorXd
          parameterVector = Map<VectorXd>(parameterArray.data(),
                                          Index(parameterArray.size()));
        }
      }
    }
  }

  return parameterVector;
}

size_t ReadParameterFromJson(
  const string &jsonFilename, 
  const string &jsonKey)
{
  // Open the JSON file
  ifstream ifs(jsonFilename);
  if (!ifs.is_open())
  {
    cerr << "Error: Unable to open file " << jsonFilename << endl;
    return 0; // Return default size_t value (0)
  }

  // Parse the JSON
  json all_parameters;
  ifs >> all_parameters;

  // Check if the main key exists
  if (all_parameters.contains(jsonKey))
  {
    const json &parameter = all_parameters[jsonKey];

    // Return the value as size_t
    return parameter.get<size_t>();
  }
  else
  {
    cerr << "Error: Key " << jsonKey << " not found in the JSON file" << endl;
    exit(3);
  }
}