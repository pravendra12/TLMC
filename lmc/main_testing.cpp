/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Pravendra                                                                       *
 * @Date: 1/16/20 3:55 AM                                                                         *
 * @Last Modified by: pravendr                                                                    *
 * @Last Modified time: 10/25/23 10:35 PM                                                         *
 **************************************************************************************************/

/*! \file  main.cpp
 *  \brief File for the main function.
 */
/*
#include "Home.h"


int main(int argc, char *argv[]) {
  if (argc == 1) {
    std::cout << "No input parameter filename." << std::endl;
    return 1;
  }
  api::Parameter parameter(argc, argv);
  api::Print(parameter);
  api::Run(parameter);
}
*/
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <iostream>

using namespace Eigen;
using namespace std;

/*
inline bool PositionCompareStateMirror111(
    const std::pair<size_t, RowVector3d> &lhs,
    const std::pair<size_t, RowVector3d> &rhs,
    const RowVector3d &mirror_direction)  // e.g., (1, 1, 1) or (-1, 1, 1)
{
  const auto &pos_lhs = lhs.second;
  const auto &pos_rhs = rhs.second;

  // Normalize the mirror direction to get the plane normal
  const RowVector3d n = mirror_direction.normalized();

  // Correct: center of the unit cube projected along the mirror direction
  const double center_proj = RowVector3d(0.5, 0.5, 0.5).dot(n);

  // Project positions onto the mirror direction (i.e., normal vector)
  const double proj_lhs = pos_lhs.dot(n) - center_proj;
  const double proj_rhs = pos_rhs.dot(n) - center_proj;

  // Primary: sort by absolute distance from the mirror plane
  const double abs_diff_proj = std::abs(proj_lhs) - std::abs(proj_rhs);
  if (abs_diff_proj < -constants::kEpsilon) { return true; }
  if (abs_diff_proj > constants::kEpsilon) { return false; }

  // Secondary: signed difference from the mirror plane
  const double signed_diff_proj = proj_lhs - proj_rhs;
  if (signed_diff_proj < -constants::kEpsilon) { return true; }
  if (signed_diff_proj > constants::kEpsilon) { return false; }

  // Tie-breaker: use x, y, z components
  const double diff_x = pos_lhs[0] - pos_rhs[0];
  if (diff_x < -constants::kEpsilon) { return true; }
  if (diff_x > constants::kEpsilon) { return false; }

  const double diff_y = pos_lhs[1] - pos_rhs[1];
  if (diff_y < -constants::kEpsilon) { return true; }
  if (diff_y > constants::kEpsilon) { return false; }

  const double diff_z = pos_lhs[2] - pos_rhs[2];
  if (diff_z < -constants::kEpsilon) { return true; }
  if (diff_z > constants::kEpsilon) { return false; }

  return false;
}
*/

/********************** Data Generation + E_KRA and Testing LCE encoder for EKRA **********************/

#include "Config.h"
#include "PotentialEnergyEstimator.h"
#include "JsonUtility.h"
#include "ClusterExpansion.h"
#include "VacancyMigrationBarrierPredictor.h"
#include "Symmetry.h"
#include "SymmetrySpglib.h"
#include "LocalEnvironmentEncoder.h"
#include "VacancyMigrationPredictor.h"

/*
RowVectorXd GetTensorProduct(const vector<RowVectorXd> &basisVector)
{
  if (basisVector.empty())
  {
    return RowVectorXd(); // return empty if no input
  }

  // Start with the first vector
  MatrixXd result = basisVector[0];

  // Iteratively compute outer products
  for (size_t i = 1; i < basisVector.size(); ++i)
  {
    result = Eigen::Map<RowVectorXd>(result.data(), result.size()).transpose() * basisVector[i];
  }

  // Flatten the final matrix row-wise
  return Eigen::Map<RowVectorXd>(result.data(), result.size());
}
*/

#include <vector>
#include <Eigen/Dense>
#include <algorithm> // for std::next_permutation

#include "AtomBasis.h"
#include "CorrelationFunction.h"
#include "GetOrbits.h"

using Eigen::RowVectorXd;
using std::vector;

template <typename T>
std::unordered_map<T, size_t, boost::hash<T>> ConvertSetToHashMapMain(
    const std::set<T> &cluster_type_set)
{
  std::unordered_map<T, size_t, boost::hash<T>> cluster_type_count;
  for (const auto &cluster_type : cluster_type_set)
  {
    cluster_type_count[cluster_type] = 0;
  }
  return cluster_type_count;
}

// returns neighbour of a given latticeId from a given sites
vector<size_t> getNeighbours(const Config &config,
                             const size_t &latticeId,
                             const vector<size_t> &allowedSites,
                             const size_t maxBondOrder)
{
  vector<size_t> neighboursLatticeIdVector;
  for (auto id : allowedSites)
  {
    auto bondOrder = config.GetDistanceOrder(latticeId, id);

    if (bondOrder > 0 && bondOrder <= maxBondOrder)
    {
      neighboursLatticeIdVector.emplace_back(id);
    }
  }

  return neighboursLatticeIdVector;
}

void getClusterCount(const Config &cfg,
                     const set<Element> elementSet,
                     const size_t maxClusterSize,
                     const size_t maxBondOrder,
                     const std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterHashSet)
{
  auto initializedClusterMap = InitializeClusterTypeSet(cfg,
                                                        elementSet,
                                                        maxClusterSize,
                                                        maxBondOrder);

  auto cluster_type_count_hashmap(ConvertSetToHashMapMain(initializedClusterMap));

  auto initializedLatticeClusterMap =
      InitializeLatticeClusterTypeSet(
          cfg,
          maxClusterSize,
          maxBondOrder);
  auto latticeClusterTypeCountMap(ConvertSetToHashMapMain(initializedLatticeClusterMap));

  for (const auto &lattice_cluster : latticeClusterHashSet)
  {
    auto atom_cluster_type = IdentifyAtomClusterType(cfg, lattice_cluster.GetLatticeIdVector());
    cluster_type_count_hashmap.at(ClusterType(atom_cluster_type, lattice_cluster.GetClusterType()))++;

    latticeClusterTypeCountMap.at(lattice_cluster.GetClusterType())++;
  }

  for (auto cluster : initializedClusterMap)
  {
    cout << cluster << " : " << cluster_type_count_hashmap.at(cluster)
         << " : " << latticeClusterTypeCountMap.at(cluster.lattice_cluster_type_) << endl;
  }
}

std::vector<std::vector<size_t>> AddOneSiteToExistingClusterHelperTest(
    const Config &reference_config,
    const size_t max_bond_order,
    const std::vector<std::vector<size_t>> &old_clusters,
    const std::vector<size_t> &allowedSites)
{

  std::vector<std::vector<size_t>> new_clusters;

  // Iterate over each existing cluster
  for (const auto &old_cluster : old_clusters)
  {
    std::set<size_t> neighbors{};

    // Retrieve neighbors up to max_bond_order
    for (auto lattice_id : old_cluster)
    {
      // allowed neighbours to define the local environment
      for (auto id : allowedSites)
      {
        auto bondOrder = reference_config.GetDistanceOrder(lattice_id, id);

        if (bondOrder > 0 && bondOrder <= max_bond_order)
        {
          neighbors.insert(id);
        }
      }
      // neighbors.insert(neighbor_lists[m][lattice_id].begin(), neighbor_lists[m][lattice_id].end());
    }

    // Remove the sites already in the cluster
    for (auto lattice_id : old_cluster)
    {
      neighbors.erase(lattice_id);
    }

    // Add new sites to the cluster
    for (auto new_lattice_id : neighbors)
    {
      std::vector<size_t> new_cluster{old_cluster};

      // Check if all the lattice sites in the old cluster are valid neighbors for the new lattice site
      bool valid_cluster = true;
      for (auto old_lattice_id : old_cluster)
      {
        auto distance_order = reference_config.GetDistanceOrder(old_lattice_id, new_lattice_id);

        // Validate distance order
        if (distance_order > max_bond_order || distance_order == 0)
        {
          valid_cluster = false;
          break;
        }
      }

      // std::cout << std::endl;

      if (valid_cluster)
      {
        new_cluster.push_back(new_lattice_id);
        new_clusters.push_back(new_cluster);
      }
    }
  }

  return new_clusters;
}

// suggestion to convert equivalentSites to vector<unordered_map>
bool isSymmetricCluster(const vector<size_t> &cluster,
                        const vector<size_t> &ssVector, // symmetrically sorted vector
                        const vector<vector<size_t>> &equivalentSites)
{
  // Build a map from site -> group_id
  // map of lattice id and in which group they lies based on 3 bar symmetry
  unordered_map<size_t, size_t> latticeIdEquivalenceMap;

  for (size_t groupId = 0; groupId < equivalentSites.size(); ++groupId)
  {
    for (auto siteId : equivalentSites[groupId])
    {
      latticeIdEquivalenceMap[ssVector[siteId]] = groupId;
    }
  }

  if (cluster.empty())
    return true;

  // Get the group of the first site
  auto groupId = latticeIdEquivalenceMap[cluster[0]];

  // Check if all sites belong to the same group
  for (auto latticeId : cluster)
  {
    if (latticeIdEquivalenceMap[latticeId] != groupId)
    {
      return false;
    }
  }

  return true;
}

std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>> FindAllEquivalentClustersTest(
    const Config &reference_config,
    size_t max_cluster_size,
    size_t max_bond_order,
    const std::pair<size_t, size_t> &latticeIdJumpPair,
    const vector<vector<size_t>> &equivalentSites)
{
  std::unordered_set<LatticeCluster, boost::hash<LatticeCluster>> lattice_cluster_hashset;

  // This is the vector which will define the local environment
  /////////////////////// need to chnage
  auto lattice_id_vector = reference_config.GetSortedLatticeVectorStateOfPair(latticeIdJumpPair, 1);

  // map for faster lookup
  // latticeId : eqSite Index
  std::unordered_map<size_t, size_t> latticeIdToIndexMap;

  size_t i = 0;
  for (auto id : lattice_id_vector)
  {
    latticeIdToIndexMap.insert(std::make_pair(id, i));
    i++;
  }

  std::vector<std::vector<size_t>> cluster_list{{}};
  if (lattice_id_vector.empty())
  {
    for (size_t i = 0; i < reference_config.GetNumLattices(); ++i)
    {
      cluster_list.push_back({i});
    }
  }
  else
  {
    for (auto lattice_id : lattice_id_vector)
    {
      cluster_list.push_back({lattice_id});
    }
  }
  int numSym = 0;
  int numNonSym = 0;

  for (size_t i = 0; i < max_cluster_size; i++)
  {
    if (i > 0)
    {
      cluster_list = AddOneSiteToExistingClusterHelperTest(reference_config, max_bond_order, cluster_list, lattice_id_vector);
    }
    for (auto &cluster : cluster_list)
    {
      auto type = IdentifyLatticeClusterType(reference_config, cluster);
      lattice_cluster_hashset.emplace(type, cluster);

      // cout << type  << endl;
      cout << "--------------------" << endl;
      if (isSymmetricCluster(cluster, lattice_id_vector, equivalentSites))
      {
        cout << "Symmetric" << endl;
        numSym++;
      }
      else
      {
        cout << "Non Symmetric" << endl;
        numNonSym++;
      }

      print1DVector(cluster);
      cout << "{ ";
      for (auto id : cluster)
      {
        cout << latticeIdToIndexMap.at(id) << " ";
      }
      cout << " }" << endl;
      cout << type << endl;

      // std::cout << type << " :  {" ;
      // for(const auto& tmp : cluster){
      //   std::cout << tmp << " " << reference_config.GetElementOfLattice(tmp) << ", ";
      // }
      // std::cout << "} : " << IdentifyAtomClusterType(reference_config,cluster) << std::endl;
    }
  }

  cout << "Total clusters: " << numSym + numNonSym << endl;
  cout << "Total symmetric clusters : " << numSym << endl;
  cout << "Total non symmetric clusters : " << numNonSym << endl;

  return lattice_cluster_hashset;
}

void GetEquivalentClusters(const Config &config,
                           const pair<size_t, size_t> &latticeJumpPair,
                           const vector<vector<size_t>> &equivalentSites)
{
  auto ssVector = config.GetSortedLatticeVectorStateOfPair(latticeJumpPair, 1);
  size_t maxBondOrder = 3;

  // Equivalent pairs

  vector<vector<size_t>> symmetricPairs;
  unordered_map<size_t, size_t> nonSymmetricPairs;

  // encoding vector for all possible pairs

  for (int i = 0; i < equivalentSites.size(); i++)
  {
    for (int j = i; j < equivalentSites.size(); j++)
    {
      if (i == j)
        cout << "Symmetric Pair" << endl;
      else
        cout << "Non Symmetric Pair" << endl;

      print1DVector(equivalentSites[i]);
      print1DVector(equivalentSites[j]);

      //

      if (i == j)
      {
        for (int m = 0; m < equivalentSites[i].size(); m++)
        {
          for (int n = m + 1; n < equivalentSites[i].size(); n++)
          {
            auto equivalentId1 = equivalentSites[i][m];
            auto equivalentId2 = equivalentSites[j][n];

            size_t id1 = ssVector[equivalentId1];
            size_t id2 = ssVector[equivalentId2];

            auto bondOrder = config.GetDistanceOrder(id1, id2);

            vector<size_t> cluster = {id1, id2};

            if (bondOrder && bondOrder <= maxBondOrder)
            {
              cout << id1 << "(" << equivalentSites[i][m] << ")" << "-"
                   << id2 << "(" << equivalentSites[i][n] << ")" << endl;

              cout << IdentifyLatticeClusterType(config, cluster) << endl;

              vector<size_t> eqClusterId = {equivalentSites[i][m], equivalentSites[i][n]};

              symmetricPairs.emplace_back(eqClusterId);
            }
          }
        }
      }
    }
  }

  print2DVector(symmetricPairs);

  // Non symmetric pairs
}
/*
void TestLCEEncodingForJump(
    RowVectorXd &encodingVectorLCE,
    const Config &cfg,
    const std::set<Element> &elementSet,
    const std::pair<int, int> &latticeJumpPair,
    const size_t maxBondOrder,
    const size_t maxClusterSize)
{
  string basisType = "Chebyshev";
  std::cout << "\n==== Testing Jump Pair: (" << latticeJumpPair.first
            << ", " << latticeJumpPair.second << ") ====\n";

  auto ssVector3Bar = GetSortedLatticeStatesForPairUnder3BarSymmetry(cfg, latticeJumpPair, maxBondOrder);
  std::cout << "Sorted Lattice States:\n";
  print1DVector(ssVector3Bar);

  auto eqSites3Bar = GetEquivalentSitesUnder3BarSymmetry(cfg, maxBondOrder, latticeJumpPair);
  std::cout << "Equivalent Sites under 3Bar Symmetry:\n";
  print2DVector(eqSites3Bar);

  auto orbitMap = GetOrbits(cfg, maxClusterSize, maxBondOrder, eqSites3Bar, ssVector3Bar);

  std::cout << "Orbit Map:\n";
  for (const auto &orbit : orbitMap)
  {
    std::cout << "----- Orbit " << orbit.first << " ------\n";
    print2DVector(orbit.second);
  }

  encodingVectorLCE = GetLocalEnvironmentEncoding(cfg, elementSet, basisType, orbitMap, ssVector3Bar);
  std::cout << "LCE Encoding:\n"
            << encodingVectorLCE << "\n";
}
            */

void processDirectorySymmetryEncodingLCE(const std::string &rootDir,
                                         const std::set<Element> &elementSet,
                                         const size_t maxBondOrder,

                                         const size_t maxBondOrderOfCluster,
                                         const size_t maxClusterSize,
                                         string &dirName)
{
  // Open output files for writing results
  const vector<double> cutoffs = {3.3, 4.7, 5.6};

  std::ofstream outFileOccupation("encodeVector_Occupation_ClusterSymmetry_ConsistentOrder_BO2" + dirName + ".txt");
  std::ofstream outFileChebyshev("encodeVector_Chebyshev_ClusterSymmetry_ConsistentOrder_BO2" + dirName + ".txt");

  if (!outFileOccupation.is_open() || !outFileChebyshev.is_open())
  {
    std::cerr << "Error opening output files for writing." << std::endl;
    return;
  }

  // Write headers for the output files
  outFileOccupation << "Folder_ID\tEncoding_Vector" << std::endl;
  outFileChebyshev << "Folder_ID\tEncoding_Vector" << std::endl;

  // Iterate over directories in rootDir
  for (const auto &entry : std::filesystem::directory_iterator(rootDir))
  {
    if (entry.is_directory())
    {
      std::string subDir = entry.path().string();
      std::string configPath = subDir + "/Config/" + dirName + "_5x5x5.cfg";

      if (std::filesystem::exists(configPath))
      {
        // Process the .cfg file
        auto cfg = Config::ReadCfg(configPath);
        cfg.UpdateNeighborList(cutoffs);

        Element vacancy("X");
        cfg.SetElementOfLattice(cfg.GetCentralAtomLatticeId(), vacancy);

        pair<size_t, size_t> latticeJumpPair = {
            cfg.GetVacancyLatticeId(),
            cfg.GetNeighborLatticeIdVectorOfLattice(cfg.GetVacancyLatticeId(), 1)[0]};

        cout << "Lattice Jump Pair: " << latticeJumpPair.first << ", " << latticeJumpPair.second << endl;

        pair<size_t, size_t> forward = latticeJumpPair;
        pair<size_t, size_t> backward = {latticeJumpPair.second,
                                         latticeJumpPair.first};

        auto eqSitesEncoding = GetEquivalentSiteEncoding3BarSymmetry(cfg, maxBondOrder);

        // auto ssVectorForward = GetSortedLatticeStatesForPairUnder3BarSymmetry(cfg, forward, maxBondOrder);
        // auto equivSitesForward = GetEquivalentSitesUnder3BarSymmetry(cfg, forward, maxBondOrder);
        auto ssVectorForward = GetSSVector3FSymmetry(cfg, forward, maxBondOrder);

        // auto ssVectorBackward = GetSortedLatticeStatesForPairUnder3BarSymmetry(cfg, backward, maxBondOrder);
        auto ssVectorBackward = GetSSVector3FSymmetry(cfg, backward, maxBondOrder);

        auto orbitMapF3Bar = GetOrbits(cfg, maxClusterSize, maxBondOrderOfCluster, ssVectorForward, eqSitesEncoding);
        auto orbitMapB3Bar = GetOrbits(cfg, maxClusterSize, maxBondOrderOfCluster, ssVectorBackward, eqSitesEncoding);

        // Occupation
        RowVectorXd encodingVectorForwardOccupation = GetLocalEnvironmentEncoding(cfg,
                                                                                  elementSet,
                                                                                  "Occupation",
                                                                                  orbitMapF3Bar);

        RowVectorXd encodingVectorBackwardOccupation = GetLocalEnvironmentEncoding(cfg,
                                                                                   elementSet,
                                                                                   "Occupation",
                                                                                   orbitMapB3Bar);

        // bool isEqualOccupation = (encodingVectorBackwardOccupation == encodingVectorForwardOccupation);
        bool isEqualOccupation = encodingVectorBackwardOccupation.isApprox(encodingVectorForwardOccupation);

        cout << "Occupation Basis Type : Are equal " << isEqualOccupation << endl;

        if (!isEqualOccupation)
        {
          break;
        }

        // Chebyshev
        RowVectorXd encodingVectorForwardChebyshev = GetLocalEnvironmentEncoding(cfg,
                                                                                 elementSet,
                                                                                 "Chebyshev",
                                                                                 orbitMapF3Bar);

        RowVectorXd encodingVectorBackwardChebyshev = GetLocalEnvironmentEncoding(cfg,
                                                                                  elementSet,
                                                                                  "Chebyshev",
                                                                                  orbitMapB3Bar);

        // bool isEqualChebyshev = (encodingVectorBackwardChebyshev == encodingVectorForwardChebyshev);
        bool isEqualChebyshev = encodingVectorBackwardChebyshev.isApprox(encodingVectorForwardChebyshev);

        cout << "Chebyshev Basis Type : Are equal " << isEqualChebyshev << endl;

        if (!isEqualChebyshev)
        {
          break;
        }

        // Write results to the respective output files
        outFileOccupation << entry.path().filename().string() << "\t";
        for (int i = 0; i < encodingVectorForwardOccupation.size(); ++i)
        {
          outFileOccupation << encodingVectorForwardOccupation[i] << " ";
        }
        outFileOccupation << std::endl;

        outFileChebyshev << entry.path().filename().string() << "\t";
        for (int i = 0; i < encodingVectorForwardChebyshev.size(); ++i)
        {
          outFileChebyshev << encodingVectorForwardChebyshev[i] << " ";
        }
        outFileChebyshev << std::endl;

        std::cout << "Processed results for " << subDir << " written to output files." << std::endl;
      }
    }
  }

  // Close the output files
  outFileOccupation.close();
  outFileChebyshev.close();
}

// PBC

void wrapMain(Vector3d &position)
{
  // Apply PBC: wrap position to [0, 1) by applying modulo operation
  position = position.array() - position.array().floor(); // Element-wise operation for all coordinates
}

void encodeVector(const vector<vector<size_t>> &eqSitesVector,
                  const vector<size_t> &ssVector)
{
  unordered_map<size_t, size_t> latticeIdIndexMap;
  for (size_t i = 0; i < ssVector.size(); i++)
  {
    latticeIdIndexMap[ssVector[i]] = i;
  }

  vector<vector<size_t>> encodingVector;

  for (const auto eqSites : eqSitesVector)
  {
    vector<size_t> subVector;
    for (const auto latticeId : eqSites)
    {
      subVector.emplace_back(latticeIdIndexMap.at(latticeId));
    }
    sort(subVector.begin(), subVector.end()); // sort each group
    encodingVector.emplace_back(subVector);
  }

  sort(encodingVector.begin(), encodingVector.end()); // sort the whole encoding vector
  print2DVector(encodingVector);
}

using namespace Eigen;

#include "KRAPredictor.h"

// Test for AtomBasis class
void testCachedAtomicBasis()
{
  // Assuming elementSet is defined (replace with your actual element set)
  set<Element> elementSet = {Element("W"), Element("Ta")}; // Example set
  string basisType = "Occupation";

  // Create AtomBasis object with element set and basis type
  AtomBasis atomicBasis(elementSet, basisType);

  // Test 1: GetCachedAtomBasis for element "W"
  // cout << "Fetching atom basis for W..." << endl;
  RowVectorXd b1 = atomicBasis.GetCachedAtomBasis(Element("W"));
  // cout << "Atom basis for W: " << b1 << endl;  // Display row vector

  // Test 2: GetCachedAtomBasis for element "Ta"
  // cout << "Fetching atom basis for Ta..." << endl;
  RowVectorXd b2 = atomicBasis.GetCachedAtomBasis(Element("Ta"));
  // cout << "Atom basis for Ta: " << b2 << endl;

  // Test 3: Call GetCachedAtomBasis for "Ta" multiple times to check caching
  // cout << "Fetching atom basis for Ta (multiple calls)..." << endl;
  atomicBasis.GetCachedAtomBasis(Element("Ta"));
  atomicBasis.GetCachedAtomBasis(Element("Ta"));
  atomicBasis.GetCachedAtomBasis(Element("Ta"));
  atomicBasis.GetCachedAtomBasis(Element("Ta"));
  atomicBasis.GetCachedAtomBasis(Element("Ta"));

  // Test 4: GetCachedAtomBasis for "W" after multiple "Ta" calls
  // cout << "Fetching atom basis for W after multiple Ta fetches..." << endl;
  atomicBasis.GetCachedAtomBasis(Element("W"));

  // Test 5: Set up basis vector and call GetCachedTensorProduct
  // cout << "Creating basis vector for tensor product..." << endl;
  vector<RowVectorXd> basisVector = {b1, b2}; // Using b1 and b2 from earlier

  string elementCluster = "WTa"; // Example cluster
  bool isSymmetric = true;       // Symmetry flag
  // cout << "Computing cached tensor product..." << endl;

  RowVectorXd tensorProduct = atomicBasis.GetCachedTensorProduct(elementCluster, basisVector, isSymmetric);
  // cout << "Computed tensor product: " << tensorProduct << endl;

  // Non symmetric
  basisVector = {b1, b2};
  isSymmetric = false;
  elementCluster = "W|Ta";

  tensorProduct = atomicBasis.GetCachedTensorProduct(elementCluster, basisVector, isSymmetric);
  // cout << "Computed tensor product: " << tensorProduct << endl;

  isSymmetric = false;
  basisVector = {b2, b1};

  elementCluster = "Ta|W";

  tensorProduct = atomicBasis.GetCachedTensorProduct(elementCluster, basisVector, isSymmetric);
  // cout << "Computed tensor product: " << tensorProduct << endl;

  isSymmetric = false;
  basisVector = {b1, b2, b2};

  elementCluster = "W|Ta|Ta";

  tensorProduct = atomicBasis.GetCachedTensorProduct(elementCluster, basisVector, isSymmetric);
  // cout << "Computed tensor product: " << tensorProduct << endl;

  // the size of tensort hash map will not change
  tensorProduct = atomicBasis.GetCachedTensorProduct(elementCluster, basisVector, isSymmetric);
  tensorProduct = atomicBasis.GetCachedTensorProduct(elementCluster, basisVector, isSymmetric);
  tensorProduct = atomicBasis.GetCachedTensorProduct(elementCluster, basisVector, isSymmetric);
}

bool CheckSymmetryEncodingConsistency(Config &cfg, vector<vector<size_t>> &expectedEncoding, size_t maxBondOrder)
{
  size_t totalSites = cfg.GetNumLattices();

  for (size_t id1 = 0; id1 < totalSites; ++id1)
  {
    auto nnSites = cfg.GetNeighborLatticeIdVectorOfLattice(id1, 1);

    for (size_t id2 : nnSites)
    {
      std::pair<size_t, size_t> forwardPair = {id1, id2};
      std::pair<size_t, size_t> backwardPair = {id2, id1};

      // ---- Forward encoding ----
      auto ssVectorF = GetSSVector3FSymmetry(cfg, forwardPair, maxBondOrder);
      // Lattice Id to index map
      std::unordered_map<size_t, size_t> mapF;
      for (size_t i = 0; i < ssVectorF.size(); ++i)
        mapF[ssVectorF[i]] = i;

      auto eqF = GetEquivalentSitesUnder3BarSymmetry(cfg, forwardPair, maxBondOrder);
      std::vector<std::vector<size_t>> encodingF;
      for (const auto &group : eqF)
      {
        std::vector<size_t> local;
        for (size_t site : group)
          local.push_back(mapF.at(site));
        std::sort(local.begin(), local.end());
        encodingF.push_back(local);
      }
      std::sort(encodingF.begin(), encodingF.end());

      // ---- Backward encoding ----
      auto ssVectorB = GetSSVector3FSymmetry(cfg, backwardPair, maxBondOrder);
      std::unordered_map<size_t, size_t> mapB;
      for (size_t i = 0; i < ssVectorB.size(); ++i)
        mapB[ssVectorB[i]] = i;

      auto eqB = GetEquivalentSitesUnder3BarSymmetry(cfg, backwardPair, maxBondOrder);
      std::vector<std::vector<size_t>> encodingB;
      for (const auto &group : eqB)
      {
        std::vector<size_t> local;
        for (size_t site : group)
          local.push_back(mapB.at(site));
        std::sort(local.begin(), local.end());
        encodingB.push_back(local);
      }
      std::sort(encodingB.begin(), encodingB.end());

      // ---- Compare ----
      if (encodingF != encodingB)
      {
        if (encodingF != expectedEncoding)
        {
          std::cout << "Mismatch found between forward and backward encodings:\n";
          std::cout << "Pair: " << id1 << " -> " << id2 << std::endl;

          std::cout << "Forward Encoding:\n";
          print2DVector(encodingF);
          std::cout << "Backward Encoding:\n";
          print2DVector(encodingB);

          std::cout << "Expected Encoding:\n";
          print2DVector(expectedEncoding);

          return false;
        }
      }
    }
  }

  std::cout << "All forward and backward encodings are consistent.\n";
  return true;
}

#include "ShortRangeOrder.h"
#include "LocalEnvironment.h"

int main()
{
  /// Implementation and validation of E_KRA model /////////////////////

  const string predictorFilename = "predictorFileKRA_BO2_WTa.json";
  const size_t maxBondOrder = 2;
  const size_t maxClusterSize = 3;
  const size_t maxBondOrderOfCluster = 3;
  /*
    /// Getting LCE
    set<Element> elementSet{Element("Ta"), Element("W")};

    vector<string> folderVector = {
        "Ta10W90",
        "Ta40W60",
        "Ta70W30",
        "Ta20W80",
        "Ta50W50",
        "Ta90W10",
        "Ta30W70",
        "Ta60W40",
    };

    for (auto folderName : folderVector)
    {
      processDirectorySymmetryEncodingLCE("/media/sf_Phd/WTaNEB/" + folderName,
                                          elementSet,
                                          maxBondOrder,
                                          maxBondOrderOfCluster,
                                          maxClusterSize,
                                          folderName);
    }

  }
  */

  const vector<double> cutoffs = {3.3, 4.7, 5.6};

  size_t increment = 1000000;
  set<Element> elementSet = {Element("Ta"), Element("W")};

  // Even though 10.7 is cutoff for the 13th NN then also it will be distance order 4
  for (int i = 0; i < 11; i++)
  {
    string filename = to_string(increment * i);

    cout << "------- " + filename + " ---------" << endl;

    auto cfg = Config::ReadCfg("//media/sf_Phd/start_50x50x50_cmc_2K_1e8.cfg.gz");

    cfg.UpdateNeighborList(cutoffs);
    

   

    
    auto start = std::chrono::high_resolution_clock::now();
    
    ShortRangeOrder sro(cfg, elementSet);
    sro.FindWarrenCowley(1);
    sro.FindWarrenCowley(2);
    sro.FindWarrenCowley(3);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Total time: " << elapsed.count() << " seconds\n";





    break;
  }

  //
  //  cout << "Central Atom Lattice ID: " << cfg.GetCentralAtomLatticeId() << endl;
  //
  //  auto nnListBO13 = cfg.GetNeighborLatticeIdVectorOfLattice(cfg.GetCentralAtomLatticeId(), 4);
  //
  //  auto nnListBO3 = cfg.GetNeighborLatticeIdVectorOfLattice(cfg.GetCentralAtomLatticeId(), 3);
  //  auto nnListBO2 = cfg.GetNeighborLatticeIdVectorOfLattice(cfg.GetCentralAtomLatticeId(), 2);
  //  auto nnListBO1 = cfg.GetNeighborLatticeIdVectorOfLattice(cfg.GetCentralAtomLatticeId(), 1);
  //
  //
  //
  //  cout << "Size of NN List 4: " << nnListBO13.size() << endl;
  //  cout << "Size of NN List 3: " << nnListBO3.size() << endl;
  //
  //
  //  // Store it
  //  auto basis = cfg.GetBasis();
  //  Eigen::Matrix3Xd relativePositionMatrix;
  //  relativePositionMatrix.resize(3, 300);
  //  vector<Element> atomVector;
  //  // atomVector.reserve(nnListBO13.size());
  //
  //  cout << "Basis" << endl;
  //  cout << basis << endl;
  //
  //  int idx = 0;
  //
  //  std::vector<size_t> allLatticeIds;
  // allLatticeIds.reserve(nnListBO13.size() + nnListBO3.size() + nnListBO2.size() + nnListBO1.size() + 1);
  //
  // allLatticeIds.insert(allLatticeIds.end(), nnListBO13.begin(), nnListBO13.end());
  // allLatticeIds.insert(allLatticeIds.end(), nnListBO3.begin(),  nnListBO3.end());
  // allLatticeIds.insert(allLatticeIds.end(), nnListBO2.begin(),  nnListBO2.end());
  // allLatticeIds.insert(allLatticeIds.end(), nnListBO1.begin(),  nnListBO1.end());
  //
  // allLatticeIds.emplace_back(cfg.GetCentralAtomLatticeId());
  //
  //
  //
  //  for (auto const latticeId : allLatticeIds)
  //  {
  //    Vector3d relativePosition = cfg.GetRelativePositionOfLattice(latticeId);
  //    auto element = cfg.GetElementOfLattice(latticeId);
  //
  //    relativePositionMatrix.col(idx) = relativePosition;
  //    atomVector.emplace_back(element);
  //    idx++;
  //  }
  //
  //  relativePositionMatrix.conservativeResize(3, idx);
  //
  //  cout << idx << endl;
  //
  //  Config supercell = Config{basis, relativePositionMatrix, atomVector};

  // Config::WriteConfig("testMFPTLocalEnv.cfg", supercell);
  // cfg.GetNeighborLatticeIdsUpToOrder(cfg.GetCentralAtomLatticeId(), 2);
  // cfg.GetNeighborLatticeIdsUpToOrder(cfg.GetCentralAtomLatticeId(), 3);
  // cfg.GetNeighborLatticeIdsUpToOrder(cfg.GetCentralAtomLatticeId(), 4);
  // cfg.GetNeighborLatticeIdsUpToOrder(cfg.GetCentralAtomLatticeId(), 5);

  /*
  auto atomVector = cfg.GetAtomVector();
  set<Element> elementSet{atomVector.begin(), atomVector.end()};

  Element vacancy("X");

  size_t vacancyId = cfg.GetCentralAtomLatticeId();
  cfg.SetElementOfLattice(vacancyId, vacancy);
  size_t nnId = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[6];

  cout << cfg.GetElementOfLattice(vacancyId).GetElementString() << " " << cfg.GetElementOfLattice(nnId).GetElementString() << endl;

  cout << vacancyId << "-" << nnId << endl;
  pair<size_t, size_t> latticeJumpPair = {vacancyId, nnId};

  KRAPredictor kraPredictor(predictorFilename, cfg, elementSet);
  kraPredictor.GetKRA(cfg, latticeJumpPair);
  kraPredictor.GetKRA(cfg, make_pair(latticeJumpPair.second, latticeJumpPair.first));

  elementSet.insert(vacancy);
  PotentialEnergyEstimator peEstimator(predictorFilename, cfg, cfg, elementSet);

  peEstimator.GetDeThreadSafe(cfg, latticeJumpPair);
  cfg.LatticeJump(latticeJumpPair);
  peEstimator.GetDeThreadSafe(cfg, latticeJumpPair);
  */

  /*
  //// Testing the LCE funciton and GetOrbits and GetEquivalentSites3BarSymmetry ///

    cout << "Testing the function" << endl;

  size_t vacancyId = 3;
  int i =0;
  auto nnSites = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1);
  for (auto id : nnSites)
  { // id = 49;


    pair<size_t, size_t> latticeJumpPair = {vacancyId,
                                            id};

    cout << latticeJumpPair.first << "-" << latticeJumpPair.second << endl;
    cout << "BO : " << cfg.GetDistanceOrder(vacancyId, id) << endl;


    // cout << "--- ---- ---- Backward ---- ---- ---- " << endl;

    // ssVector = GetSortedVector3BarSymmetry(cfg, make_pair(latticeJumpPair.second, latticeJumpPair.first), 1);

    // print1DVector(ssVector);


    cout << "Symmetrically sorted under 3 Bar symmetry" << endl;

    auto ssVector = GetSSVector3FSymmetry(cfg,latticeJumpPair, 2);

    print1DVector(ssVector);

    unordered_map<size_t, size_t> latticeToIndexMap;


    size_t localIdx = 0;
    for (size_t localId : ssVector)
    {
      latticeToIndexMap[localId] = localIdx;
      localIdx++;
    }


    cout << "----- Equivalent Sites -----" << endl;
    auto eqSites3Bar = GetEquivalentSitesUnder3BarSymmetry(cfg, latticeJumpPair, 2);

    print2DVector(eqSites3Bar);

    cout << " ----- Encoded ----- " <<endl;
    vector<vector<size_t>> encoding;

    for (auto eqSite : eqSites3Bar)
    {
      vector<size_t> localEncoding;
      for (auto site : eqSite)
      {
        localEncoding.emplace_back(latticeToIndexMap.at(site));
      }
      sort(localEncoding.begin(), localEncoding.end());
      encoding.emplace_back(localEncoding);
    }

    print2DVector(encoding);

    cout << "--- ---- ---- Backward ---- ---- ---- " << endl;
    cout << "Symmetrically sorted under 3 Bar symmetry" << endl;

    ssVector = GetSSVector3FSymmetry(cfg,make_pair(latticeJumpPair.second, latticeJumpPair.first), 2);

    print1DVector(ssVector);

    size_t idx = 0;
    for (size_t id : ssVector)
    {
      latticeToIndexMap[id] = idx;
      idx++;
    }


    cout << "----- Equivalent Sites -----" << endl;
    eqSites3Bar = GetEquivalentSitesUnder3BarSymmetry(cfg, make_pair(latticeJumpPair.second, latticeJumpPair.first), 2);

    print2DVector(eqSites3Bar);

    cout << " ----- Encoded ----- " <<endl;
    vector<vector<size_t>> encodingB;

    for (auto eqSite : eqSites3Bar)
    {
      vector<size_t> localEncoding;
      for (auto site : eqSite)
      {
        localEncoding.emplace_back(latticeToIndexMap.at(site));
      }
      sort(localEncoding.begin(), localEncoding.end());
      encodingB.emplace_back(localEncoding);
    }

    print2DVector(encodingB);


     //
     // cout << "----- Orbit Map -----" << endl;
     // for (const auto &orbit : orbitMap)
     // {
     //   std::cout << "----- Orbit " << orbit.first << " ------\n";
     //   // print2DVector(orbit.second);
     //   vector<vector<size_t>> orbitEncoding;
     //   for (auto latticeIdVector : orbit.second)
     //   {
     //    vector<size_t> localEncoding;
     //    for (auto latticeId : latticeIdVector)
     //    {
     //      localEncoding.emplace_back(latticeToIndexMap.at(latticeId));
     //    }
     //    sort(localEncoding.begin(), localEncoding.end());
     //    orbitEncoding.emplace_back(localEncoding);
     //   }
     //   print2DVector(orbit.second);
     //   cout << " ----- Orbit Encoding ----- " << endl;
     //   print2DVector(orbitEncoding);
     // }



    }

    auto eqSiteEncoding = GetEquivalentSiteEncoding3BarSymmetry(cfg, 2);

    print2DVector(eqSiteEncoding);

    /// Checking whether the encodings are consistent for all the possible jump pairs
    /*
    vector<vector<size_t>> expectedEncodingBO2 = {{0, 19, },
    {7, 8, 9, 10, 11, 12, },
    {4, 5, 6, 13, 14, 15, },
    {1, 2, 3, 16, 17, 18, }};



    vector<vector<size_t>> expectedEncodingBO3  {{0, 1, 2, 35, 36, 37, },
      {3, 34, },
      {4, 5, 6, 31, 32, 33, },
      {7, 8, 9, 28, 29, 30, },
      {10, 13, 14, 23, 25, 27, },
      {11, 12, 15, 22, 24, 26, },
      {16, 17, 18, 19, 20, 21, }};


    // For BO = 2
    CheckSymmetryEncodingConsistency(cfg, expectedEncodingBO3, 3);
    */
}