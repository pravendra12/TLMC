#include "JsonUtility.h"

pair<VectorXd, double>
ReadParametersFromJson(const string &json_filename,
                       const string &json_key)
{
  ifstream ifs(json_filename, ifstream::in);
  json all_parameters;
  ifs >> all_parameters;

  VectorXd adjustedBetaVector;
  double adjustedIntercept;

  for (const auto &[key, parameters] : all_parameters.items())
  {
    if (key == json_key)
    {
      for (const auto &[subKey, subParams] : parameters.items())
      {
        if (subKey == "adjusted_beta_" + json_key)
        {

          vector<double> adjustedBetaVec;

          for (const auto &adjustedBeta : subParams)
          {
            adjustedBetaVec.push_back(adjustedBeta.get<double>());
          }

          // Convert the vector to VectorXd
          adjustedBetaVector = Map<VectorXd>(adjustedBetaVec.data(),
                                             Index(adjustedBetaVec.size()));
        }
        else if (subKey == "adjusted_intercept_" + json_key)
        {
          adjustedIntercept = subParams[0].get<double>();
        }
      }
    }
  }

  return make_pair(move(adjustedBetaVector), adjustedIntercept);
}

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
