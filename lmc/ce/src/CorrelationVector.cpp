#include "CorrelationVector.h"
#include "PrintUtility.h"

// Make use of symmetry
/*
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
*/
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

// GetCorrelationVector with a element at specific site
VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const size_t &targetLatticeId,
    const Element &elementToAssign,
    const vector<size_t> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<size_t>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<size_t> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        targetLatticeId,
        elementToAssign,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

VectorXd GetCorrelationVector(
    const Config &config,
    BasisSet &atomicBasis,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentOrbitVector)
{
  VectorXd correlationVector = VectorXd::Zero(equivalentOrbitVector.size());

  for (size_t idx = 0; idx < equivalentOrbitVector.size(); ++idx)
  {
    const auto &orbitPair = equivalentOrbitVector[idx];

    // Handle empty cluster
    if (orbitPair.first.size() == 1 && orbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        config,
        atomicBasis,
        orbitPair.first);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

double GetOrbitCorrelationFunction(const Config &config,
                                   BasisSet &atomicBasis,
                                   const vector<vector<size_t>> &orbitVector)
{

  double orbitClusterFunction = 0;

  for (const auto &latticeCluster : orbitVector)
  {
    vector<Element> elementVector;
    elementVector.reserve(latticeCluster.size());

    for (const auto &latticeId : latticeCluster)
    {
      auto element = config.GetElementOfLattice(latticeId);
      elementVector.emplace_back(element);
    }

    AtomClusterType atomClusterType(elementVector);

    double clusterFunction = atomicBasis.GetBasisProduct(atomClusterType);
    orbitClusterFunction += clusterFunction;
  }

  orbitClusterFunction = orbitClusterFunction / static_cast<double>(orbitVector.size());

  return orbitClusterFunction;
}

// Assiging element at a given targetLatticeId
double GetOrbitCorrelationFunction(const Config &config,
                                   BasisSet &atomicBasis,
                                   const size_t &targetLatticeId,
                                   const Element &elementToAssign,
                                   const vector<vector<size_t>> &orbitVector)
{

  double orbitClusterFunction = 0;

  for (const auto &latticeCluster : orbitVector)
  {
    vector<Element> elementVector;
    elementVector.reserve(latticeCluster.size());

    for (const auto &latticeId : latticeCluster)
    {
      // Assigning the element to a specific site
      if (latticeId == targetLatticeId)
      {
        elementVector.emplace_back(elementToAssign);
      }
      else
      {
        auto element = config.GetElementOfLattice(latticeId);
        elementVector.emplace_back(element);
      }
    }

    AtomClusterType atomClusterType(elementVector);

    double clusterFunction = atomicBasis.GetBasisProduct(atomClusterType);
    orbitClusterFunction += clusterFunction;
  }

  orbitClusterFunction = orbitClusterFunction / static_cast<double>(orbitVector.size());

  return orbitClusterFunction;
}

/// FOR TLMC

double GetOrbitCorrelationFunction(
    const TiledSupercell &tiledSupercell,
    const size_t &smallConfigIdx, // smallConfigIdx
    BasisSet &atomicBasis,
    const vector<size_t> &neighborsOfSmallConfig, // neighbours of those smallConfig
    const vector<vector<LatticeSiteEncodedMapping>> &orbitVector)
{

  double orbitClusterFunction = 0;

  for (const auto &latticeCluster : orbitVector)
  {
    vector<Element> elementVector;
    elementVector.reserve(latticeCluster.size());

    for (const auto &latticeSiteEncodedMapping : latticeCluster)
    {
      auto siteLatticeId = latticeSiteEncodedMapping.latticeId;
      int encodedSmallConfigIdx = latticeSiteEncodedMapping.encodedSmallConfigId;

      size_t mappedSmallConfigIdx;
      if (encodedSmallConfigIdx == -1)
      {
        mappedSmallConfigIdx = smallConfigIdx; // same small config index or current config Index
      }
      else
      {
        mappedSmallConfigIdx = neighborsOfSmallConfig[encodedSmallConfigIdx];
      }
      // Create a lattice-site mapping
      LatticeSiteMapping siteMapping{siteLatticeId, mappedSmallConfigIdx};

      // Get the element at that site
      const auto &element = tiledSupercell.GetElementAtSite(siteMapping);
      elementVector.emplace_back(element);
    }

    AtomClusterType atomClusterType(elementVector);

    double clusterFunction = atomicBasis.GetBasisProduct(atomClusterType);
    orbitClusterFunction += clusterFunction;
  }

  orbitClusterFunction = orbitClusterFunction / static_cast<double>(orbitVector.size());

  return orbitClusterFunction;
}

VectorXd GetCorrelationVector(
    const TiledSupercell &tiledSupercell,
    const LatticeSiteMapping &latticeSite, // Lattice Site around which correlation vector is computed
    BasisSet &atomicBasis,
    const vector<LatticeSiteEncodedMapping> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  const Cube &cubeObj = tiledSupercell.GetCube();
  // Neighbours of the smallConfig in which latticeId lies
  const vector<size_t> &neighborsOfSmallConfig = cubeObj.GetNeighbors(latticeSite.smallConfigId);

  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<LatticeSiteEncodedMapping>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<LatticeSiteEncodedMapping> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        tiledSupercell,
        latticeSite.smallConfigId,
        atomicBasis,
        neighborsOfSmallConfig,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

/// FOR ELEMENT AT SPECIFIC SITE

// Assiging element at a given targetLatticeId
double GetOrbitCorrelationFunction(
    const TiledSupercell &tiledSupercell,
    const size_t &smallConfigIdx, // smallConfigIdx
    BasisSet &atomicBasis,
    const LatticeSiteMapping &targetLatticeSite,
    const Element &elementToAssign,
    const vector<size_t> &neighborsOfSmallConfig, // neighbours of those smallConfig
    const vector<vector<LatticeSiteEncodedMapping>> &orbitVector)
{

  double orbitClusterFunction = 0;

  for (const auto &latticeCluster : orbitVector)
  {
    vector<Element> elementVector;
    elementVector.reserve(latticeCluster.size());

    for (const auto &latticeSiteEncodedMapping : latticeCluster)
    {
      auto siteLatticeId = latticeSiteEncodedMapping.latticeId;
      int encodedSmallConfigIdx = latticeSiteEncodedMapping.encodedSmallConfigId;

      size_t mappedSmallConfigIdx;
      if (encodedSmallConfigIdx == -1)
      {
        mappedSmallConfigIdx = smallConfigIdx; // same small config index or current config Index
      }
      else
      {
        mappedSmallConfigIdx = neighborsOfSmallConfig[encodedSmallConfigIdx];
      }
      // Create a lattice-site mapping
      LatticeSiteMapping siteMapping{siteLatticeId, mappedSmallConfigIdx};

      // Assigning the element to a specific site
      if (siteMapping == targetLatticeSite)
      {
        elementVector.emplace_back(elementToAssign);
      }
      else
      {
        const auto &element = tiledSupercell.GetElementAtSite(siteMapping);
        elementVector.emplace_back(element);
      }
    }

    AtomClusterType atomClusterType(elementVector);

    double clusterFunction = atomicBasis.GetBasisProduct(atomClusterType);
    orbitClusterFunction += clusterFunction;
  }

  orbitClusterFunction = orbitClusterFunction / static_cast<double>(orbitVector.size());

  return orbitClusterFunction;
}

// GetCorrelationVector with a element at specific site
VectorXd GetCorrelationVector(
    const TiledSupercell &tiledSupercell,
    const LatticeSiteMapping &latticeSite, // Lattice Site around which correlation vector is computed
    BasisSet &atomicBasis,
    const LatticeSiteMapping &targetLatticeSite,
    const Element &elementToAssign,
    const vector<LatticeSiteEncodedMapping> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{
  const Cube &cubeObj = tiledSupercell.GetCube();
  // Neighbours of the smallConfig in which latticeId lies
  const vector<size_t> &neighborsOfSmallConfig = cubeObj.GetNeighbors(latticeSite.smallConfigId);

  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<LatticeSiteEncodedMapping>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<LatticeSiteEncodedMapping> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        tiledSupercell,
        latticeSite.smallConfigId,
        atomicBasis,
        targetLatticeSite,
        elementToAssign,
        neighborsOfSmallConfig,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

//// Correlation Vector For KRA Predictor
// Need to think about better designing this class

VectorXd GetCorrelationVector(
    const TiledSupercell &tiledSupercell,
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair,
    BasisSet &atomicBasis,
    const vector<NeighbourOfPair> &canonicalSortedLatticeIds,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &equivalentEncodedClusters)
{

  VectorXd correlationVector = VectorXd::Zero(equivalentEncodedClusters.size());

  for (size_t idx = 0; idx < equivalentEncodedClusters.size(); ++idx)
  {
    const auto &encodedOrbitPair = equivalentEncodedClusters[idx];

    // Handle empty cluster
    if (encodedOrbitPair.first.size() == 1 && encodedOrbitPair.first[0].empty())
    {
      correlationVector(idx) = 1.0; // φ₀ = 1
      continue;
    }

    vector<vector<NeighbourOfPair>> orbitVector;
    orbitVector.reserve(encodedOrbitPair.first.size());

    for (const auto &encodedCluster : encodedOrbitPair.first)
    {
      vector<NeighbourOfPair> cluster;
      cluster.reserve(encodedCluster.size());

      for (auto encodedIdx : encodedCluster)
      {
        cluster.emplace_back(
            canonicalSortedLatticeIds[encodedIdx]);
      }
      orbitVector.emplace_back(cluster);
    }

    double orbitCorrelationFunction = GetOrbitCorrelationFunction(
        tiledSupercell,
        latticeSiteJumpPair,
        atomicBasis,
        orbitVector);

    correlationVector(idx) = orbitCorrelationFunction;

    // cout << orbitPair.second << " : " << orbitCorrelationFunction << endl;
  }

  return correlationVector;
}

double GetOrbitCorrelationFunction(
    const TiledSupercell &tiledSupercell,
    const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteJumpPair,
    BasisSet &atomicBasis,
    const vector<vector<NeighbourOfPair>> &orbitVector)
{
  const Cube &cubeObj = tiledSupercell.GetCube();

  // Neighbours of the smallConfig in which latticeId lies
  // For Site 1
  const vector<size_t> &neighborsOfSmallConfig1 = cubeObj.GetNeighbors(latticeSiteJumpPair.first.smallConfigId);

  // For Site 2
  const vector<size_t> &neighborsOfSmallConfig2 = cubeObj.GetNeighbors(latticeSiteJumpPair.second.smallConfigId);

  double orbitClusterFunction = 0;

  for (const auto &latticeCluster : orbitVector)
  {
    vector<Element> elementVector;
    elementVector.reserve(latticeCluster.size());

    for (const auto &neighbourOfPairEncoded : latticeCluster)
    {
      const auto &siteMappingInfo = neighbourOfPairEncoded.latticeSiteInfo;
      const size_t siteLatticeId = siteMappingInfo.latticeId;
      const int encodedSmallConfigIdx = siteMappingInfo.encodedSmallConfigId;
      // If this neighbourOfPairEncoded belongs to first latticeId then use
      // neighbours of first else second also smallConfigIdx

      LatticeSiteMapping siteMapping;
      if (neighbourOfPairEncoded.origin == NeighbourOfPair::SourceSite::First)
      {
        size_t mappedSmallConfigIdx;
        if (encodedSmallConfigIdx == -1) // use the base site itself
        {
          mappedSmallConfigIdx = latticeSiteJumpPair.first.smallConfigId; // same small config index or current config Index
        }
        else
        {
          mappedSmallConfigIdx = neighborsOfSmallConfig1[encodedSmallConfigIdx];
        }
        // Create a lattice-site mapping
        siteMapping.latticeId = siteLatticeId;
        siteMapping.smallConfigId = mappedSmallConfigIdx;
      }

      else if (neighbourOfPairEncoded.origin == NeighbourOfPair::SourceSite::Second)
      {
        size_t mappedSmallConfigIdx;
        if (encodedSmallConfigIdx == -1) // use the base site itself
        {
          mappedSmallConfigIdx = latticeSiteJumpPair.second.smallConfigId; // same small config index or current config Index
        }
        else
        {
          mappedSmallConfigIdx = neighborsOfSmallConfig2[encodedSmallConfigIdx];
        }
        // Create a lattice-site mapping
        siteMapping.latticeId = siteLatticeId;
        siteMapping.smallConfigId = mappedSmallConfigIdx;
      }

      // Get the element at that site
      const auto &element = tiledSupercell.GetElementAtSite(siteMapping);
      elementVector.emplace_back(element);
    }

    AtomClusterType atomClusterType(elementVector);

    double clusterFunction = atomicBasis.GetBasisProduct(atomClusterType);
    orbitClusterFunction += clusterFunction;
  }

  orbitClusterFunction = orbitClusterFunction / static_cast<double>(orbitVector.size());

  return orbitClusterFunction;
}
