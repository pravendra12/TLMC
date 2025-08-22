/*******************************************************************************
 * Copyright (c) 2025. All rights reserved.
 * @Author: Pravendra Patel
 * @Date:    2025-06-01
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 *******************************************************************************/

/*! @file KRAPredictor.cpp
    @brief File contains implementation of Kinetically Resolved Activation Barrier
 */

#include "KRAPredictor.h"

KRAPredictor::KRAPredictor(
    const ClusterExpansionParameters &ceParams,
    const Config &config) : maxBondOrder_(ceParams.GetMaxBondOrderKRA()),
                            maxBondOrderOfCluster_(ceParams.GetMaxBondOrderOfClusterKRA()),
                            maxClusterSize_(ceParams.GetMaxClusterSizeKRA()),
                            atomicBasis_(
                                ceParams.GetElementSetKRA(),
                                ceParams.GetBasisType()),
                            kecisMap_(ceParams.GetKECIs()),
                            canonicalReferenceMap_(
                                GetCenteredNeighborsAlongJumpDirection(
                                    config,
                                    maxBondOrder_,
                                    referenceJumpDirection_)),
                            symmetryOperations_(GetSymmetryOperations(config)),
                            equivalentClustersEncoding_(
                                GetEquivalentClustersEncoding(
                                    config,
                                    maxBondOrder_,
                                    maxBondOrderOfCluster_,
                                    maxClusterSize_,
                                    canonicalReferenceMap_))

{
}

double KRAPredictor::GetKRA(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair)
{
  const auto &elementFirst = config.GetElementOfLattice(latticeIdJumpPair.first);
  const auto &elementSecond = config.GetElementOfLattice(latticeIdJumpPair.second);

  // At least one element must be a vacancy
  if (elementFirst != Element("X") && elementSecond != Element("X"))
  {
    throw runtime_error("Error in `KRAPredictor::GetKRA`: Jump pair must contain at least one vacancy.");
  }

  Element migratingElement = (elementFirst == Element("X"))
                                 ? elementSecond
                                 : elementFirst;

  // Canonical Pair
  size_t id1 = latticeIdJumpPair.first;
  size_t id2 = latticeIdJumpPair.second;

  pair<size_t, size_t> canonicalJumpPair = (id1 < id2)
                                               ? std::make_pair(id1, id2)
                                               : std::make_pair(id2, id1);

  auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForPair(
      config,
      canonicalJumpPair,
      maxBondOrder_,
      canonicalReferenceMap_,
      symmetryOperations_);

  VectorXd correlationVector = GetCorrelationVector(
      config,
      atomicBasis_,
      canonicalSortedLatticeIds,
      equivalentClustersEncoding_);

  // E = J.Φ_α

  VectorXd kecis = kecisMap_.at(migratingElement);
  double eKraValue = kecis.dot(correlationVector);

  return eKraValue;
}
