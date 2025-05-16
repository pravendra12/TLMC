#include "LocalEnvironmentEncoder.h"

/**
 * @brief Checks if the given orbit ID represents a symmetric orbit.
 *
 * This function determines whether all characters in the provided orbit ID
 * string are identical. An empty orbit ID is considered symmetric by default.
 *
 * @param orbitId The orbit ID string to check.
 * @return true if the orbit ID is empty or all characters in the string are the same;
 *         false otherwise.
 */
inline bool isSymmetricOrbit(const string &orbitId)
{
  if (orbitId.empty())
    return true;

  return all_of(orbitId.begin(), orbitId.end(), [&](char c)
                { return c == orbitId[0]; });
}

// This function concatenate vector which contanin RowVectorXd in to RowVectorXd
static RowVectorXd ConcatenateVector(const std::vector<RowVectorXd> &vectors)
{
  std::vector<double> buffer;
  buffer.reserve(std::accumulate(vectors.begin(), vectors.end(), static_cast<size_t>(0),
                                 [](size_t sum, const RowVectorXd &v)
                                 { return sum + v.size(); }));

  for (const auto &v : vectors)
  {
    buffer.insert(buffer.end(), v.data(), v.data() + v.size());
  }

  RowVectorXd result = Eigen::Map<RowVectorXd>(buffer.data(), buffer.size());
  return result;
}

RowVectorXd GetLocalEnvironmentEncoding(
    const Config &config,
    const set<Element> &elementSet,
    const string &basisType,
    const map<string, vector<vector<size_t>>> &orbitEncodingMap)
{
  vector<RowVectorXd> encodingVector;
  encodingVector.reserve(orbitEncodingMap.size());

  // Iterate over the orbitEncodingMap
  for (const auto &orbitEncoding : orbitEncodingMap)
  {
    bool isClusterSymmetric = isSymmetricOrbit(orbitEncoding.first);

    // Get Correlation function for each orbit
    RowVectorXd orbitCorrelationFunction = GetCorrelationFunction(config,
                                                                  elementSet,
                                                                  basisType,
                                                                  orbitEncoding.second,
                                                                  isClusterSymmetric);

    if (orbitCorrelationFunction.size() > 0)
      encodingVector.emplace_back(orbitCorrelationFunction);
  }

  return ConcatenateVector(encodingVector);
}