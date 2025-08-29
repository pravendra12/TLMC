#include "CorrelationVector.h"

// Make use of symmetry
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  vector<VectorXd> correlationFunctions;
  correlationFunctions.reserve(equivalentEncodedClusters.size());

  for (const auto &encodedOrbitPair : equivalentEncodedClusters)
  {
    
    if (encodedOrbitPair.first.size() == 1)
    {
      if (encodedOrbitPair.first[0].empty())
      {
        VectorXd orbitCorrelationFunction(1);
        orbitCorrelationFunction(0) = 1.0; // φ₀ = 1

        correlationFunctions.emplace_back(orbitCorrelationFunction);

        continue;
      }
    }

    VectorXd orbitCorrelationFunction = GetCorrelationFunction(
        config,
        atomicBasis,
        canonicalSortedLatticeIds,
        encodedOrbitPair.first);

    
    if (orbitCorrelationFunction.size() > 0)
      correlationFunctions.emplace_back(orbitCorrelationFunction);
  }

  // Concat the correlationFunctions and return a VectorXd

  size_t totalSize = 0;
  for (const auto &v : correlationFunctions)
    totalSize += v.size();

  VectorXd correlationVector(totalSize);

  size_t offset = 0;
  for (const auto &v : correlationFunctions)
  {
    correlationVector.segment(offset, v.size()) = v;
    offset += v.size();
  }

  return correlationVector;
}