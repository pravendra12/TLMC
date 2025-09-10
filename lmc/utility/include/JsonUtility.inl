#include "JsonUtility.h"

template <typename T>
T ReadParameterFromJson(
    const json &parameters,
    const string &jsonKey)
{
  // Check if the main key exists
  if (parameters.contains(jsonKey))
  {
    const json &parameter = parameters[jsonKey];
    return parameter.get<T>();
  }
  else
  {
    cerr << "Error: JSON key '" << jsonKey << "' not found in the coefficient file" << endl;
    throw runtime_error("JSON key '" + jsonKey + "' not found");
  }
}


template <typename T>
T ReadSingleParameterFromJson(const json &allParameters,
                              const string &jsonKey,
                              const string &jsonSubKey)
{
    if (!allParameters.contains(jsonKey))
        throw std::runtime_error("JSON key not found");

    const json &parameters = allParameters[jsonKey];

    if (!parameters.contains(jsonSubKey))
        throw std::runtime_error("JSON sub-key not found");

    return parameters[jsonSubKey].get<T>();
}


template <typename T>
vector<T> ReadParametersFromJson(const json &allParameters,
                                 const string &jsonKey,
                                 const string &jsonSubKey)
{
  vector<T> parameterVector;

  // Check if the main key exists
  if (!allParameters.contains(jsonKey))
  {
    cerr << "Error: JSON key '" << jsonKey << "' not found in the coefficient file" << endl;
    throw runtime_error("JSON key '" + jsonKey + "' not found");
  }

  const json &parameters = allParameters[jsonKey];

  // Check if the sub key exists
  if (!parameters.contains(jsonSubKey))
  {
    cerr << "Error: JSON sub-key '" << jsonSubKey << "' not found under '" << jsonKey << "'" << endl;
    throw runtime_error("JSON sub-key '" + jsonSubKey + "' not found");
  }

  const json &subParams = parameters[jsonSubKey];

  // Extract values
  for (const auto &param : subParams)
  {
    parameterVector.push_back(param.get<T>());
  }

  return parameterVector;
}
