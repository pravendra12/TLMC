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

void processDirectorySymmetryEncodingLCE(const std::string &rootDir,
                                         const std::set<Element> &elementSet,
                                         const size_t maxBondOrder,

                                         const size_t maxBondOrderOfCluster,
                                         const size_t maxClusterSize,
                                         string &dirName)
{
  // Open output files for writing results
  const vector<double> cutoffs = {3.3, 4.7, 5.6};

  std::ofstream outFileOccupation("encodeVector_Occupation_ClusterSymmetry" + dirName + ".txt");
  std::ofstream outFileChebyshev("encodeVector_Chebyshev_ClusterSymmetry" + dirName + ".txt");

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

        auto ssVectorForward = GetSortedLatticeStatesForPairUnder3BarSymmetry(cfg, forward, maxBondOrder);
        auto equivSitesForward = GetEquivalentSitesUnder3BarSymmetry(cfg, maxBondOrder, forward);

        auto ssVectorBackward = GetSortedLatticeStatesForPairUnder3BarSymmetry(cfg, backward, maxBondOrder);
        auto equivSitesBackward = GetEquivalentSitesUnder3BarSymmetry(cfg, maxBondOrder, backward);

        auto orbitMapF3Bar = GetOrbits(cfg, maxClusterSize, maxBondOrderOfCluster, equivSitesForward, ssVectorForward);
        auto orbitMapB3Bar = GetOrbits(cfg, maxClusterSize, maxBondOrderOfCluster, equivSitesBackward, ssVectorBackward);

        // Occupation
        RowVectorXd encodingVectorForwardOccupation = GetLocalEnvironmentEncoding(cfg,
                                                                                  elementSet,
                                                                                  "Occupation",
                                                                                  orbitMapF3Bar,
                                                                                  ssVectorForward);

        RowVectorXd encodingVectorBackwardOccupation = GetLocalEnvironmentEncoding(cfg,
                                                                                   elementSet,
                                                                                   "Occupation",
                                                                                   orbitMapB3Bar,
                                                                                   ssVectorBackward);

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
                                                                                 orbitMapF3Bar,
                                                                                 ssVectorForward);

        RowVectorXd encodingVectorBackwardChebyshev = GetLocalEnvironmentEncoding(cfg,
                                                                                  elementSet,
                                                                                  "Chebyshev",
                                                                                  orbitMapB3Bar,
                                                                                  ssVectorBackward);

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



using namespace Eigen;

#include "KRAPredictor.h"

int main()
{
  /// Implementation and validation of E_KRA model /////////////////////


  const string predictorFilename = "predictorFileEKRA_WTa.json";
  const size_t maxBondOrder = 3;
  const size_t maxClusterSize = 3;
  const size_t maxBondOrderOfCluster = 3;

  const vector<double> cutoffs = {3.3, 4.7, 5.6};

  auto cfg = Config::ReadCfg("/media/sf_Phd/WTaNEB/Ta90W10/10/Config/Ta90W10_5x5x5.cfg");
  cfg.UpdateNeighborList(cutoffs);



  Element vacancy("X");

  size_t vacancyId = cfg.GetCentralAtomLatticeId();

  cfg.SetElementOfLattice(vacancyId, vacancy);


  auto atomVector = cfg.GetAtomVector();
  set<Element> elementSet{atomVector.begin(), atomVector.end()};

  pair<size_t, size_t> latticeJumpPair = {vacancyId,
                                          cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[0]};

  cout << latticeJumpPair.first << endl;
  cout << latticeJumpPair.second << endl;


  // KRAPredictor kraPredictor(predictorFilename,cfg,  elementSet, maxBondOrder, maxClusterSize);

  // kraPredictor.GetKRA(cfg, latticeJumpPair);

  // kraPredictor.GetKRA(cfg, make_pair(latticeJumpPair.second, latticeJumpPair.first));

  // elementSet.insert(Element("X"));

  // 
  // PotentialEnergyEstimator peEstimator(predictorFilename,
  //                                      cfg,
  //                                      cfg,
  //                                      elementSet,
  //                                      maxClusterSize,
  //                                      maxBondOrder);
// 
  // double dE = peEstimator.GetDeThreadSafe(cfg, latticeJumpPair);
  // cout << dE << endl;
  
  // cfg.LatticeJump(latticeJumpPair);
  // dE = peEstimator.GetDeThreadSafe(cfg, make_pair(latticeJumpPair.second, latticeJumpPair.first));
  // cout << dE << endl;


  VacancyMigrationPredictor predictor(predictorFilename, 
                                      cfg,
                                      cfg, 
                                      elementSet,maxClusterSize, maxBondOrder);


  auto dEAndBarrier = predictor.GetBarrierAndDeltaE(cfg, latticeJumpPair);
  cout <<"Barrier : " << dEAndBarrier.first << endl;
  cout <<"dE : " << dEAndBarrier.second << endl;

  cfg.LatticeJump(latticeJumpPair);

  dEAndBarrier = predictor.GetBarrierAndDeltaE(cfg, latticeJumpPair);
  cout <<"Barrier : " << dEAndBarrier.first << endl;
  cout <<"dE : " << dEAndBarrier.second << endl;

  cout << cfg.GetVacancyLatticeId() << endl;






  







  




  ///////////////// Generate E_KRA Data ////////////////////////
  /*
  // Process the directory
  size_t maxClusterSize = 3;
  size_t maxBondOrder = 3;
  size_t maxBondOrderOfCluster = 3;

  set<Element> elementSet{Element("Ta"), Element("W")};

  // processDirectory("/media/sf_Phd/WTaNEB/nebTa50W50/Ta50W50/Ta50W50/",
  //                  elementSet,
  //                  maxBondOrderOfCluster,
  //                  maxClusterSizeOfCluster);

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
  */
  /*
  string predictorFilename = "predictor_file_WTa_version2.json";
  const vector<double> cutoffs = {3.3, 4.7, 5.6};

  auto cfg = Config::ReadCfg("/media/sf_Phd/WTaNEB/Ta90W10/05/Config/Ta90W10_5x5x5.cfg");
  cfg.UpdateNeighborList(cutoffs);

  Element vacancy("X");

  size_t vacancyId = cfg.GetCentralAtomLatticeId();

  cfg.SetElementOfLattice(vacancyId, vacancy);

  pair<size_t, size_t> latticeJumpPair = {vacancyId,
                                          cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[0]};

  cout << latticeJumpPair.first << endl;
  cout << latticeJumpPair.second << endl;

  auto atomVector = cfg.GetAtomVector();
  set<Element> elementSet{atomVector.begin(), atomVector.end()};

  RowVectorXd encodingLCEForward;
  RowVectorXd encodingLCEBackward;

  size_t maxBondOrder = 3;
  size_t maxClusterSize = 3;

  pair<size_t, size_t> latticeJumpPair1 = make_pair(13, cfg.GetNeighborLatticeIdVectorOfLattice(13, 1)[0]);

  // cout << cfg.GetRelativePositionOfLattice()

  // auto eqSites = GetEquivalentSitesUnder3BarSymmetry(cfg, 1, latticeJumpPair);
  // print2DVector(eqSites);
  //
  // eqSites = GetEquivalentSitesUnder3BarSymmetry(cfg, 1, make_pair(latticeJumpPair.second,
  //                                                                 latticeJumpPair.first));
  // print2DVector(eqSites);
  //
  // eqSites = GetEquivalentSitesUnder3BarSymmetry(cfg, 1, latticeJumpPair1);
  // print2DVector(eqSites);
  //
  // eqSites = GetEquivalentSitesUnder3BarSymmetry(cfg, 1, make_pair(latticeJumpPair1.second, latticeJumpPair1.first));
  // print2DVector(eqSites);

  // auto orbitMap = GetOrbits(cfg, 2, 1, )

  cout << "==== Central Pair ======" << endl;
  cout << "==== Forward Jump ======" << endl;
  TestLCEEncodingForJump(encodingLCEForward,
                         cfg,
                         elementSet,
                         latticeJumpPair,
                         maxBondOrder,    // maxBondOrder
                         maxClusterSize); // maxClusterSize

  cout << "==== Backward Jump ======" << endl;
  TestLCEEncodingForJump(encodingLCEBackward,
                         cfg,
                         elementSet,
                         make_pair(latticeJumpPair.second, latticeJumpPair.first),
                         maxBondOrder,    // maxBondOrder
                         maxClusterSize); // maxClusterSize

  bool isEqual = (encodingLCEBackward == encodingLCEForward);
  cout << "Are equal : " << isEqual << endl;

  cout << "==== Boundary Pair ======" << endl;

  cout << "==== Forward Jump ======" << endl;
  TestLCEEncodingForJump(encodingLCEForward,
                         cfg,
                         elementSet,
                         latticeJumpPair1,
                         maxBondOrder,    // maxBondOrder
                         maxClusterSize); // maxClusterSize

  cout << "==== Backward Jump ======" << endl;
  TestLCEEncodingForJump(encodingLCEBackward,
                         cfg,
                         elementSet,
                         make_pair(latticeJumpPair1.second, latticeJumpPair1.first),
                         maxBondOrder,    // maxBondOrder
                         maxClusterSize); // maxClusterSize
  cout << "Are equal : " << isEqual << endl;
  */
  

  // auto nnFirst = cfg.GetNeighboringLatticeIdSetOfPair(latticeJumpPair, 1);

  // for (auto id : nnFirst)
  // {
  //   cout << id << " : " << cfg.GetRelativePositionOfLattice(id).transpose() << endl;
  // }

  // Testing the new symmetry fuction

  /*
  auto supercellCfg = Config::GenerateSupercell(5, 3.243, "X", "BCC");
  supercellCfg.UpdateNeighborList(cutoffs);

  PotentialEnergyEstimator peEstimator(predictorFilename,
                                       cfg,
                                       supercellCfg,
                                       elementSet, 2, 1);

  auto equivalentSites = GetEquivalentSitesUnderKFoldRotation(cfg, 1, 3);

  print2DVector(equivalentSites);

  auto ssVector = cfg.GetSortedLatticeVectorStateOfPair(latticeJumpPair, 1);
  cout << "{ ";
  for (auto id : ssVector)
  {
    cout << id << cfg.GetElementOfLattice(id) << ", ";
  }
  cout << " }" << endl;

  // print1DVector(ssVector);

  vector<size_t> cluster = {latticeJumpPair.first, latticeJumpPair.second};

  // peEstimator.GetEncodeVectorOfCluster(cfg, cluster);

  auto allLatticeClusters = FindAllLatticeClusters(cfg, 2, 1, cluster);

  for (auto latticeCluster : allLatticeClusters)
  {
    cout << "---------------------------------\n"
         << endl;
    print1DVector(latticeCluster.GetLatticeIdVector());
    cout << IdentifyAtomClusterType(cfg, latticeCluster.GetLatticeIdVector()) << endl;
    cout << latticeCluster.GetClusterType() << endl;
  }

  vector<size_t> cluster1 = {latticeJumpPair.first, latticeJumpPair.second};
  cout << "Cluster Forward: " << endl;
  peEstimator.GetEncodeVectorOfCluster(cfg, cluster1);

  cfg.LatticeJump(latticeJumpPair);

  // vector<size_t> cluster2 = {latticeJumpPair.second, latticeJumpPair.first};
  cout << "Cluster Backward: " << endl;
  peEstimator.GetEncodeVectorOfCluster(cfg, cluster1);

  // Symmetry operation

  // RowVector3d pos(0, 0, 0);
  // auto dirs = Get111Directions();
  //
  // // cout << "Normalized"  << pos.normalized() << endl;
  //
  // for (const auto& dir : dirs) {
  //   RowVector3d mirror = MirrorAcrossDirection(pos.normalized(), dir);
  //   std::cout << "Mirror across " << dir
  //                 << " -> " << mirror.transpose() << "\n";
  // }

  cout << cfg.GetRelativePositionOfLattice(latticeJumpPair.first).transpose() << endl;
  cout << cfg.GetRelativePositionOfLattice(latticeJumpPair.second).transpose() << endl;

  // direction

  RowVector3d vec1 = cfg.GetRelativePositionOfLattice(latticeJumpPair.first).transpose();
  RowVector3d vec2 = cfg.GetRelativePositionOfLattice(latticeJumpPair.second).transpose();

  RowVector3d n = (vec1 + vec2) / 2;
  RowVector3d n_caps = n.normalized();

  cout << n << endl;

  for (auto id : ssVector)
  {
    cout << id << " " << cfg.GetRelativePositionOfLattice(id).transpose() << endl;
  }

  auto eqSites = GetEquivalentSitesUnderKFoldRotation(cfg, 1, 3);

  print2DVector(eqSites);

  vector<size_t> allCluster{ssVector.begin(), ssVector.end()};
  allCluster.emplace_back(cluster[0]);
  allCluster.emplace_back(cluster[1]);

  auto allLatticeClustersOverall = FindAllLatticeClusters(cfg, 2, 1, allCluster);

  int i = 0;
  for (auto latticeCluster : allLatticeClustersOverall)
  {
    cout << "---------------------------------\n"
         << endl;
    print1DVector(latticeCluster.GetLatticeIdVector());
    cout << IdentifyAtomClusterType(cfg, latticeCluster.GetLatticeIdVector()) << endl;
    cout << latticeCluster.GetClusterType() << endl;

    i += 1;
  }
  cout << "Total clusters: " << i << endl;

  // Under 3 bar symmetry

  vector<vector<size_t>> eqSites3Bar = {{0, 13}, {1, 2, 3, 11, 12, 10}, {4, 5, 6, 7, 8, 9}};

  // singlets easy
  /*
  cfg.LatticeJump(latticeJumpPair);

  for (int i = 0; i < eqSites3Bar.size(); i++)
  {
    // Iterate over
    for (int j = i; j < eqSites3Bar.size(); j++)
    {
      auto eqSitesVector1 = eqSites3Bar[i];
      auto eqSitesVector2 = eqSites3Bar[j];

      print1DVector(eqSitesVector1);
      print1DVector(eqSitesVector2);

      // Try to make pair between them
      cout << "--------------------\n"
           << endl;

      for (int id1 = 0; id1 < eqSitesVector1.size(); id1++)
      {
        for (int id2 = 0; id2 < eqSitesVector2.size(); id2++)
        {
          if (cfg.GetDistanceOrder(allCluster[eqSitesVector1[id1]], allCluster[eqSitesVector2[id2]]) < 3 &&
              cfg.GetDistanceOrder(allCluster[eqSitesVector1[id1]], allCluster[eqSitesVector2[id2]]))
          {
            vector<size_t> vec{allCluster[eqSitesVector1[id1]], allCluster[eqSitesVector2[id2]]};
            cout << eqSitesVector1[id1] << "-" << eqSitesVector2[id2] << " "
                 << IdentifyLatticeClusterType(cfg, vec) << " " << IdentifyAtomClusterType(cfg, vec)
                 << " " << cfg.GetElementOfLattice(vec[0]) << " " << cfg.GetElementOfLattice(vec[1]) << endl;
          }
        }
      }

      cout << "--------------------\n"
           << endl;
    }
  }
  // cout << "Total clusters: " << i << endl;

  */
  /*
  GetEquivalentClusters(cfg, latticeJumpPair, eqSites3Bar);
  */
  // Using existing funciton to get the type of cluster then change the IdentifyAtomClusterType
  // based on the atom arrangemetn and symmetry but for latticeClusterType we can define
  // latticeClsuterType map for a given local enviroment

  // Not a good idea
  // FindAllEquivalentClusters(cfg, 2, 1, latticeJumpPair);

  // testing add one site to existing cluster
  // vector<vector<size_t>> clusterSmall = {{132, 136}};
  //
  // auto newCluster = AddOneSiteToExistingClusterHelperTest(cfg,
  //                                                         3, clusterSmall, ssVector);
  // print2DVector(newCluster);
  //
  // for (auto cluster : newCluster)
  // {
  //   cout << IdentifyLatticeClusterType(cfg, cluster) << endl;
  // }
  // Works well  and good
  /*
  // will return pairs
  elementSet.erase(Element("X"));
  auto latticeClusterHashSet = FindAllEquivalentClustersTest(cfg, 3, 3, latticeJumpPair, eqSites3Bar);

  getClusterCount(cfg, elementSet,
                  3, 3, latticeClusterHashSet);

  auto ssVectorF = cfg.GetSortedLatticeVectorStateOfPair(latticeJumpPair, 1);
  auto ssVectorB = cfg.GetSortedLatticeVectorStateOfPair(
      make_pair(latticeJumpPair.second,
                latticeJumpPair.first),
      1);

  cout << "-------------------------" << endl;
  cout << "Forward" << endl;
  getClusterCount(cfg, elementSet, 3, 3, FindClustersWithinAllowedSites(cfg, 3, 3, ssVectorF));
  */
  /*
  cout << "-------------------------" << endl;
  cout << "Backward" << endl;
  getClusterCount(cfg, elementSet, 3, 3, FindClustersWithinAllowedSites(cfg, 3, 3, ssVectorB));
  */

  /*
  const Config &reference_config,
  const size_t max_bond_order,
  const std::vector<std::vector<size_t>> &old_clusters,
  const std::vector<size_t> &allowedSites
  */
  /*
  unordered_map<size_t, size_t> latticeIdToIndexMap;
  for (size_t i = 0; i < ssVectorF.size(); i++)
  {
    latticeIdToIndexMap.insert(make_pair(ssVectorF[i], i));
  }

  unordered_map<size_t, size_t> indexToOrbitMap;

  for (size_t i = 0; i < eqSites3Bar.size(); i++)
  {
    for (auto j : eqSites3Bar[i])
    {
      indexToOrbitMap.insert(make_pair(j, i));
    }
  }

  auto firstNNLatticeCluster = FindClustersWithinAllowedSites(cfg, 2, 1, ssVectorF);

  unordered_map<string, vector<vector<size_t>>> orbitMap;

  for (auto cluster : firstNNLatticeCluster)
  {
    auto latticeIdVector = cluster.GetLatticeIdVector();

    // get the index and its orbit

    cout << "----------------------" << endl;
    print1DVector(latticeIdVector);
    vector<size_t> indexVector;
    vector<size_t> orbitVector;

    // corresopind index vector and orbit

    int orbitSum = 0;
    for (auto id : latticeIdVector)
    {
      auto index = latticeIdToIndexMap.at(id);
      auto orbit = indexToOrbitMap.at(index);

      indexVector.emplace_back(index);
      orbitVector.emplace_back(orbit);

      orbitSum += orbit;
    }
    print1DVector(indexVector);
    print1DVector(orbitVector);

    // Sort orbits to create a unique key
    sort(orbitVector.begin(), orbitVector.end());
    string orbitKey;
    for (auto o : orbitVector)
    {
      orbitKey += to_string(o);
    }

    cout << orbitKey << endl;

    // Handle empty clusters separately
    if (indexVector.empty())
    {
      orbitMap["-1"].emplace_back(indexVector);
    }
    else
    {
      orbitMap[orbitKey].emplace_back(indexVector);
    }
  }

  cout << "**********" << endl;

  vector<vector<size_t>> orbit12;

  for (auto orbitInfo : orbitMap)
  {
    cout << "---------------------" << endl;
    cout << " Orbit : " << orbitInfo.first << endl;

    // compute correlation function

    print2DVector(orbitInfo.second);

    if (orbitInfo.first == "222")
    {
      orbit12 = orbitInfo.second;
    }
  }

  getClusterCount(cfg, elementSet, 2, 1, firstNNLatticeCluster);

  // apply_symmetry(cfg, requiredCluster);

  // apply_symmetry_cluster(cfg, requiredCluster);

  auto basis1 = GetAtomBasis(Element("W"), elementSet, "Occupation");
  auto basis2 = GetAtomBasis(Element("Ta"), elementSet, "Occupation");

  // Tensor product

  // W-W
  cout << "W-W" << endl;
  cout << basis1 << endl;
  cout << GetTensorProduct(basis1, basis1) << endl;

  // Ta-Ta

  cout << "Ta-Ta" << endl;
  cout << basis2 << endl;
  cout << GetTensorProduct(basis2, basis2) << endl;

  // Ta-W

  cout << "Ta-W" << endl;
  cout << GetTensorProduct(basis2, basis1) << endl;

  // W-Ta

  cout << "W-Ta" << endl;
  cout << GetTensorProduct(basis1, basis2) << endl;

  string basisType = "Chebyshev";

  // orbit 12;
  // correlation fuction
  RowVectorXd corrFunction;
  bool isCorrFunctionResized = false;
  int numClusters = 0;
  for (auto latticeCluster : orbit12)
  {
    vector<RowVectorXd> atomBasisVector;

    cout << "-------------------" << endl;
    string elementCluster = "";
    for (auto idx : latticeCluster)
    {
      auto element = cfg.GetElementOfLattice(ssVectorF[idx]);

      elementCluster += element.GetElementString();

      RowVectorXd atomBasis = GetAtomBasis(element,
                                           elementSet,
                                           basisType);

      atomBasisVector.emplace_back(atomBasis);

      cout << element.GetElementString() << " " << atomBasis << endl;
    }

    cout << elementCluster << endl;

    RowVectorXd clusterBasisVector = GetTensorProduct(atomBasisVector, true);

    cout << clusterBasisVector << endl;

    if (!isCorrFunctionResized)
    {
      corrFunction.resize(clusterBasisVector.size());
      corrFunction.setZero();
      isCorrFunctionResized = true;
    }

    corrFunction += clusterBasisVector;
    numClusters++;
  }

  cout << corrFunction / numClusters << endl;

  cout << "From the fuction" << endl;
  GetCorrelationFunction(cfg,
                         elementSet,
                         basisType,
                         orbit12,
                         true,
                         ssVectorB);


  // Testing GetLocalEnvironmentEncoder Function
  cout << "// Testing GetLocalEnvironementEncoder Function" << endl;
  RowVectorXd encodingBackward = GetLocalEnvironmentEncoding(cfg,
  elementSet,
  basisType,
  orbitMap,
  ssVectorB);


  RowVectorXd encodingForward = GetLocalEnvironmentEncoding(cfg,
    elementSet,
    basisType,
    orbitMap,
    ssVectorF);



  auto eqSites3BarTest = GetEquivalentSitesUnderKBarSymmetry(cfg, 1, 3);

  print2DVector(eqSites3BarTest);

  auto newOrbitMap = GetOrbits(cfg, 2, 1, eqSites3Bar, ssVectorF);

  for (auto orbitInfo : newOrbitMap)
  {
    cout << "---------------------" << endl;
    cout << " Orbit : " << orbitInfo.first << endl;

    // compute correlation function

    print2DVector(orbitInfo.second);
  }


  cout << "Backward Encoding : " << encodingBackward << endl;
  cout << "Forward Encoding : " << encodingForward << endl;
  */
}

/*
/********************** Testing barrier and ce version 2 **********************/
/*
#include "Config.h"
#include "PotentialEnergyEstimator.h"
#include "JsonUtility.h"
#include "ClusterExpansion.h"
#include "VacancyMigrationBarrierPredictor.h"
#include "Symmetry.h"

using namespace Eigen;
int main()
{
  string predictorFilename = "predictor_file_WTa_version2.json";

  // // Eigen::VectorXd betaCE = ReadParametersFromJson(predictorFilename,
  //                                                 "ce",
  //                                                 "beta_ce");
//
  // cout << "Beta CE: " << betaCE << endl;
//
  // VectorXd betaBarrier = ReadParametersFromJson(predictorFilename,
  //                                               "barrier",
  //                                               "beta_barrier");
//
  // cout << "Beta Barrier: " << betaBarrier << endl;
//
  // VectorXd meanXFeatures = ReadParametersFromJson(predictorFilename,
  //                                                 "barrier",
  //                                                 "mean_X_features");
//
  // VectorXd stdXFeatures = ReadParametersFromJson(predictorFilename,
  //                                                "barrier",
  //                                                "std_X_features");
//
  // VectorXd meanYBarrier = ReadParametersFromJson(predictorFilename,
  //                                                "barrier",
  //                                                "mean_Y_barrier");
//
  // VectorXd stdYBarrier = ReadParametersFromJson(predictorFilename,
  //                                               "barrier",
  //                                               "std_Y_barrier");
//
  // cout << "Mean Y Barrier: " << meanYBarrier << endl;
  // cout << "Std Y Barrier: " << stdYBarrier << endl;
//
  // // Testing and predicting the values for a know configuration and how close it
  // // is to the value obtained from LAMMPS
//
  auto cfg = Config::ReadConfig("testingVersion2/01 Ta50W50_5x5x5.cfg");

  vector<double> cutoffs = {3.22, 4.5, 5.3};

  cfg.UpdateNeighborList(cutoffs);

  auto supercellCfg = Config::GenerateSupercell(5, 3.24, "X", "BCC");
  supercellCfg.UpdateNeighborList(cutoffs);

  // {vacancyId, migratingAtomId}

  // dE estimation

  // Initial configuration
  // cout << cfg.GetCartesianPositionOfLattice(137) << endl;

  Element vacancy("X");
  cfg.SetElementOfLattice(137, vacancy);

  // Previous Version
  // dE computation

  auto atomVector = cfg.GetAtomVector();
  set<Element> elementSet{atomVector.begin(), atomVector.end()};

  PotentialEnergyEstimator peEstimator(predictorFilename,
                                       cfg,
                                       supercellCfg,
                                       elementSet,
                                       3, 3);

  // cout << "dE: \n" << peEstimator.GetDeThreadSafe(cfg, latticeIdJumpPair) << endl;


  // The prevous CE seems to be more accurate as it was trained on a variety of
  // data points

  /// Barrier Prediction


  VacancyMigrationBarrierPredictor barrierPredictor(cfg, elementSet, predictorFilename);

  // Upto 2nd NN
  // this gives backward barrier that in prevously the migration atom need to be first
  // followed by the vacancy so this is the backward barrier
  // cout << "Barrier Version 1 Backward: " << barrierPredictor.GetBarrier(cfg, latticeIdJumpPair) << endl;

  // cout << "Barrier Version 1 Forward: " << barrierPredictor.GetBarrier(cfg, latticeIdJumpPairF) << endl;
  //
  // cout << endl;




  // New barrier
  // Improved accuracy
  // Upto 3rd NN
  pair<size_t, size_t> latticePair = {cfg.GetVacancyLatticeId(), 112};
  cfg.LatticeJump(latticePair);
  for (auto id : cfg.GetNeighborLatticeIdVectorOfLattice(cfg.GetVacancyLatticeId(),1))
  {

    pair<size_t, size_t> latticeIdJumpPair = {cfg.GetVacancyLatticeId(), id};
    pair<size_t, size_t> latticeIdJumpPairBackward = {id, cfg.GetVacancyLatticeId()};

    cout << cfg.GetVacancyLatticeId()  << " -> " << id << endl;



    cout << cfg.GetElementOfLattice(latticeIdJumpPair.first).GetElementString() << " "
        << cfg.GetElementOfLattice(latticeIdJumpPair.second).GetElementString() << endl;

    cout << "Forward Barrier: " << barrierPredictor.GetBarrier(cfg, latticeIdJumpPair) << endl;


    cout << cfg.GetElementOfLattice(latticeIdJumpPairBackward.first).GetElementString() << " "
        << cfg.GetElementOfLattice(latticeIdJumpPairBackward.second).GetElementString() << endl;
    cout << "Backward Barrier: " << barrierPredictor.GetBarrier(cfg, latticeIdJumpPairBackward) << endl;

    cout << "dE: " << peEstimator.GetDeThreadSafe(cfg, latticeIdJumpPair) << endl;

    cout << "--------------------------------------------------------------" << endl;
  }

  // In simulation


  Config::WriteConfig("start_Ta50W50_5x5x5.cfg.gz", cfg);

  }

  // Result of the test

  // Previous Version

  /*
  01
  dE = -0.197

  Barrier
  Forward : 2.06615
  Backward : 2.01324

  */

// Updated Version

/*
dE = -0.2352

Forward = 2.16667
Backward = 2.431
*/
/*
// Summary

So, new version of the barrier predictor is much better that the previous version
as it is only trained on the 50-50 composition but the dE does not perform well
enough but previous dE estimator is okay, but caution need to be taken as
adjusted beta is being used but it would be more better to scale the encodings.

// Need to work on the dE part

}

/**************************** Computing Energy of Config **********************/
/*
#include "Config.h"
#include "VacancyMigrationBarrierPredictor.h"
#include <chrono>
#include "PotentialEnergyEstimator.h"
#include "EncodingUtility.h"
#include "PrintUtility.h"
#include "VacancyMigrationPredictor.h"
#include "ClusterDynamics.h"

int main(int argc, char* argv[])
{
  // Testing

  if (argc == 1) {
    std::cout << "No input config file. Usage lmc.exe <configFile>" << std::endl;
    return 1;
  }

  string filename = argv[1];

  auto cfg = Config::ReadCfg(filename);

  cfg.UpdateNeighborList({3.3, 4.7, 5.6});

  size_t maxBondOrder = 3;
  size_t maxClusterSize = 3;

  auto supercellCfg = Config::GenerateSupercell(5, 3.4, "X", "BCC");
  supercellCfg.UpdateNeighborList({3.3, 4.7, 5.6});

  auto atomVector = cfg.GetAtomVector();
  set<Element> elementSet(atomVector.begin(), atomVector.end());

  PotentialEnergyEstimator peEstimator("predictor_file_WTa_version2.json",
                                       cfg,
                                       supercellCfg,
                                       elementSet,
                                       3, 3);

  auto energyOfConfig = peEstimator.GetEnergy(cfg);

  cout << "Energy of the configuration: " << energyOfConfig << endl;


}

/*******************************ClusterDynamicsTesting*************************/
/*
#include "Config.h"
#include "VacancyMigrationBarrierPredictor.h"
#include <chrono>
#include "PotentialEnergyEstimator.h"
#include "EncodingUtility.h"
#include "PrintUtility.h"
#include "VacancyMigrationPredictor.h"
#include "ClusterDynamics.h"

int main()
{
  // Testing

  auto cfg = Config::ReadCfg("b2Cluster/10000000_1200K.cfg");

  cfg.UpdateNeighborList({3.3, 4.7, 5.6});

  // clusterSize(cfg)
  ClusterDynamics clusterDynamics(cfg);

  std::map<std::string, Config::VectorVariant> auxiliaryList;

  vector<int> clusterSizeVector;

  clusterDynamics.detectB2Clusters(auxiliaryList, clusterSizeVector);

  map<string, Config::ValueVariant> globalList;

  globalList["b2OrderParam"] = 0.45;

  cout << "Cluster Type from main: " << endl;

  std::vector<int> clusterTypes = std::get<std::vector<int>>(auxiliaryList["clusterType"]);
  int i = 0;
  for (int id : clusterTypes)
  {
    std::cout << i << " : " << clusterTypes[i] << " " << id << std::endl;
    i++;
  }

  // Cluster size

  for (auto clusterSize : clusterSizeVector)
  {
    cout << clusterSize << ", ";
  }
  cout << endl;

  // Config::WriteXyzExtended("b2Cluster/10000000_1200K.xyz.gz", cfg, auxiliaryList, globalList);
}

/*
int main()
{
  auto cfg = Config::ReadCfg("start_W50Ta50_20x20x20.cfg");

  cfg.UpdateNeighborList({3.3, 4.7, 5.6});

  auto supercellCfg = Config::GenerateSupercell(5, 3.4, "X", "BCC");

  supercellCfg.UpdateNeighborList({3.3, 4.7, 5.6});


  auto atomVector = cfg.GetAtomVector();
  set<Element> elementSet(atomVector.begin(), atomVector.end());

  VacancyMigrationBarrierPredictor barrier_predictor(cfg,
                                                     elementSet,
                                                     "predictor_file_WTa.json");

  PotentialEnergyEstimator peEstimator("predictor_file_WTa.json",
                                       cfg,
                                       supercellCfg,
                                       elementSet,
                                       3, 3);

  size_t vacancyId = cfg.GetVacancyLatticeId();
  size_t nnAtomId = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[0];

  pair<size_t, size_t> jumpPairForward = {vacancyId, nnAtomId};

  auto dE = peEstimator.GetDeThreadSafe(cfg, jumpPairForward);

  cout << dE << endl;


  pair<size_t, size_t> jumpPairBackward = {nnAtomId, vacancyId};

  auto dE_b = peEstimator.GetDeThreadSafe(cfg, jumpPairBackward);

  cout << dE_b << endl;


  // VancanMigrationPredictor

  VacancyMigrationPredictor vacancyMigrationPredictor("predictor_file_WTa.json",
                                       cfg,
                                       supercellCfg,
                                       elementSet,
                                       3, 3);

  //

  auto barrierDeForward =   vacancyMigrationPredictor.getBarrierAndEnergyChange(cfg,
                                                                         jumpPairForward);


auto barrierDeBackward =   vacancyMigrationPredictor.getBarrierAndEnergyChange(cfg,
                                                                         jumpPairBackward);

  cout << "Forward" << endl;

  cout << barrierDeForward.first << endl;
  cout << barrierDeForward.second << endl;


  cout << "Backward" << endl;

  cout << barrierDeBackward.first << endl;
  cout << barrierDeBackward.second << endl;










}

// // #include "ClusterExpansion.h"
// // #include "PotentialEnergyEstimator.h"
// using namespace std;
//
// // // #include <chrono>
//
// static std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>> ConvertLatticeSetToHashMap(
//     const std::set<LatticeClusterType> &lattice_cluster_type_set) {
//   std::unordered_map<LatticeClusterType, size_t, boost::hash<LatticeClusterType>> lattice_cluster_type_count;
//   for (const auto &cluster_type : lattice_cluster_type_set) {
//     lattice_cluster_type_count[cluster_type] = 0;
//   }
//   return lattice_cluster_type_count;
// }
//
// void PrintNeighborLists(const Config& config) {
//     const auto& neighborLists = config.GetNeighborLists();
//
//     for (size_t i = 0; i < neighborLists.size(); ++i) {
//         std::cout << "Layer " << i << ":\n";
//         for (size_t j = 0; j < neighborLists[i].size(); ++j) {
//             std::cout << "  Row " << j << ": ";
//             for (size_t k = 0; k < neighborLists[i][j].size(); ++k) {
//                 std::cout << neighborLists[i][j][k] << " ";
//             }
//             std::cout << "\n";
//         }
//     }
// }
//
// void PrintAtomVector(const Config& config) {
//     const auto& atomVector = config.GetAtomVector();
//
//     std::cout << "Atom Vector:\n";
//     for (const auto& atom : atomVector) {
//         std::cout << atom << " ";  // Adjust to `atom.GetString()` if `Element` doesn't support direct printing.
//     }
//     std::cout << "\n";
// }

// // #include "ThermodynamicAveraging.h"
// // #include "CanonicalMcSerial.h"
// // #include "CanonicalMcAbstract.h"
// // #include "PotentialEnergyEstimator.h"
// // #include "Config.h"
// // #include "ClusterExpansion.h"
// // #include <unordered_set>
// // #include <unordered_map>
// // #include <vector>
// // #include <string>
// // #include <set>
// // #include <algorithm>
// // #include <boost/functional/hash.hpp>
// // #include "CanonicalMcSerial.h"
// // #include "Symmetry.h"
// // #include <omp.h>
// // #include <chrono>
// // #include <Eigen/Dense>
// // #include <vector>
// // #include <algorithm>
// // // #include "KineticMcChainOmpi.h"
// // // #include "KineticMcAbstract.h"
// // #include "JumpEvent.h"
// // #include "Constants.hpp"
// // #include "ShortRangeOrder.h"
// // #include "Traverse.h"
// // #include "KineticMcFirstMpi.h"
// // #include "LatticeClusterMMM.hpp"
// // #include "Element.hpp"
// // #include <cmath>

//
// std::vector<size_t> GetSymmetricallySortedLatticeVectorMMM(
//     const Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair, std::vector<size_t> lattice_id_hashset) {
// Number of first-, second-, and third-nearest neighbors of the jump pairs
// constexpr size_t kNumOfSites = constants::kNumThirdNearestSetSizeOfPair;

// Get the set of neighbor lattice IDs
// auto lattice_id_hashset = config.GetNeighborsLatticeIdSetOfPair(lattice_id_jump_pair);

// Calculate the movement vector to center the pair
//  Eigen::Vector3d move_distance = Eigen::Vector3d(0.5, 0.5, 0.5);
//     // - config.GetLatticePairCenter(lattice_id_jump_pair);
//
//  // Prepare the lattice list
//  //std::vector<cfg::Lattice> lattice_list;
//  //lattice_list.reserve(kNumOfSites);
//
//  for (const auto id : lattice_id_hashset) {
//
//    auto relative_position = config.GetRelativePositionOfLattice(id);
//
//    // Move to the center
//    relative_position += move_distance;
//
//    // Ensure relative position stays within bounds (0, 1) using Eigen's .array() operations
//    relative_position = relative_position.array() - relative_position.array().floor();
//
//    std::cout << relative_position << std::endl;
//
//    //lattice.SetRelativePosition(relative_position);
//
//    }
//
//    //lattice_list.push_back(lattice);
//  }

//   // Apply rotation matrix to the lattice vector
//   Eigen::Matrix3d rotation_matrix = config.GetLatticePairRotationMatrix(lattice_id_jump_pair);
//   for (auto &lattice : lattice_list) {
//     Eigen::Vector3d relative_position = lattice.GetRelativePosition();
//     relative_position = rotation_matrix * relative_position;
//     relative_position = relative_position.array() - relative_position.array().floor();
//     lattice.SetRelativePosition(relative_position);
//   }
//
//   // Sort the lattice list based on the custom comparison function
//   std::sort(lattice_list.begin(), lattice_list.end(), [](const cfg::Lattice &lhs, const cfg::Lattice &rhs) -> bool {
//     return PositionCompareMMM(lhs, rhs);
//   });
//
//   return lattice_list;
// }

// std::unordered_set<size_t> getUniqueNeighbors(Config &cfg, size_t lattice_id1, size_t lattice_id2) {
//     std::unordered_set<size_t> unique_neighbors;
//
//         // Add neighbors for lattice_id1 up to the 3rd nearest neighbor, excluding lattice_id1 and lattice_id2
//     for (int nn_order = 1; nn_order <= 3; ++nn_order) {
//         auto neighbors = cfg.GetNeighborLatticeIdVectorOfLattice(lattice_id1, nn_order);
//         for (const auto& neighbor : neighbors) {
//             if (neighbor != lattice_id1 && neighbor != lattice_id2) {
//                 unique_neighbors.insert(neighbor);
//             }
//         }
//     }
//
//     // Add neighbors for lattice_id2 up to the 3rd nearest neighbor, excluding lattice_id1 and lattice_id2
//     for (int nn_order = 1; nn_order <= 3; ++nn_order) {
//         auto neighbors = cfg.GetNeighborLatticeIdVectorOfLattice(lattice_id2, nn_order);
//         for (const auto& neighbor : neighbors) {
//             if (neighbor != lattice_id1 && neighbor != lattice_id2) {
//                 unique_neighbors.insert(neighbor);
//             }
//         }
//     }
//
//     return unique_neighbors;
// }

// // #include <iostream>
// // #include "Config.h"
// // #include "PotentialEnergyEstimator.h"
// // #include "VacancyMigrationBarrierPredictor.h"
// // #include "JsonUtility.h"
// // #include "Symmetry.h"
// // #include "PrintUtility.h"
// // #include "Home.h"
// using namespace std;

// #include "Config.h"
// #include "Home.h"
// #include "Symmetry.h"
// #include <boost/functional/hash.hpp>
// #include <chrono>
// #include <omp.h>
// #include <mutex>
// #include <unordered_set>
// #include "PrintUtility.h"
// #include "B2OrderParameter.h"

/*
int main(int argc, char* argv[])
{
  const vector<double> cutoffs = {3.3, 4.7, 5.6};
  auto cfg = Config::ReadPoscar("NbTaW_5x5x5.POSCAR");
  cfg.UpdateNeighborList(cutoffs);


  Element ele1("Cl");
  Element ele2("Cs");

  Element eleW("W");
  Element eleTa("Ta");


  B2OrderParameter b2OrderParamWTa(cfg);

  cout << "W at Alpha : " << b2OrderParamWTa.GetAlphaSiteOccupancy(eleW) << endl;
  cout << "Ta at Alpha : " << b2OrderParamWTa.GetAlphaSiteOccupancy(eleTa) << endl;


  cout << "W at beta : " << b2OrderParamWTa.GetBetaSiteOccupancy(eleW) << endl;
  cout << "Ta at beta : " << b2OrderParamWTa.GetBetaSiteOccupancy(eleTa) << endl;


  cout << "b2 Order param for W : " << b2OrderParamWTa.GetB2OrderParameter(eleW) << endl;
  cout << "b2 Order param for Ta : " << b2OrderParamWTa.GetB2OrderParameter(eleTa) << endl;


  auto aSites = b2OrderParamWTa.GetAlphaLatticeSites();
  auto bSites = b2OrderParamWTa.GetBetaLatticeSites();

  for (auto id : aSites)
  {
    cfg.SetElementOfLattice(id, ele1);
  }

  for (auto id : bSites)
  {
    cfg.SetElementOfLattice(id, ele2);
  }

  Config::WriteConfig("CsCl_B2_Test.cfg", cfg);

  B2OrderParameter b2OrderCsCl(cfg);

  cout << "ele1 at Alpha : " << b2OrderCsCl.GetAlphaSiteOccupancy(ele1) << endl;
  cout << "ele2 at Alpha : " << b2OrderCsCl.GetAlphaSiteOccupancy(ele2) << endl;


  cout << "ele1 at beta : " << b2OrderCsCl.GetBetaSiteOccupancy(ele1) << endl;
  cout << "ele2 at beta : " << b2OrderCsCl.GetBetaSiteOccupancy(ele2) << endl;

  cout << "b2 ele1 : " << b2OrderCsCl.GetB2OrderParameter(ele1) << endl;
  cout << "b2 ele2 : " << b2OrderCsCl.GetB2OrderParameter(ele2) << endl;















}
*/
/*
int main()
{

  const vector<double> cutoffs = {3.3, 4.7, 5.6};
  auto cfg = Config::ReadCfg("start_W50Ta50_20x20x20.cfg");
  cfg.UpdateNeighborList(cutoffs);

  bool isSame = cfg.GetNumAtoms() == cfg.GetNumLattices();

  std::cout << isSame << std::endl;


  unordered_map<pair<size_t, size_t>,
  vector<size_t>,
  boost::hash<std::pair<size_t, size_t>>> symmetricallySortedVectorMap_;

  auto start = chrono::high_resolution_clock::now();

  for (size_t id1 = 0; id1 < cfg.GetNumLattices(); ++id1)
  {
    auto id1FirstNN = cfg.GetNeighborLatticeIdVectorOfLattice(id1, 1);
    for (auto &id2 : id1FirstNN)
    {
      std::pair<size_t, size_t> latticePair = {id1, id2};
      auto symmetricallySortedVector = GetSortedLatticeVectorStateOfPair(cfg, latticePair, 2);
      symmetricallySortedVectorMap_[latticePair] = symmetricallySortedVector;

    }
  }

  auto end = chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsed = end - start;

  std::cout << "Total time required : " << elapsed << std::endl;


  unordered_map<pair<size_t, size_t>,
  vector<size_t>,
  boost::hash<std::pair<size_t, size_t>>> symmetricallySortedVectorMap2;

  std::cout << "Number of threads : " << omp_get_num_threads() << std::endl;

 auto start1 = std::chrono::high_resolution_clock::now();

 ComputeSymmetricallySortedVectorMap(cfg, 3, symmetricallySortedVectorMap2);

 auto end1 = std::chrono::high_resolution_clock::now();

 std::chrono::duration<double> elapsed1 = end1 - start1;

 std::cout << "Time Taken : " << elapsed1 << std::endl;



}

/*
int main()
{
  const vector<double> cutoffs = {3.3, 4.7, 5.6};
  auto cfg = Config::ReadCfg("start_W50Ta50_20x20x20.cfg");
  cfg.UpdateNeighborList(cutoffs);

  auto supercellCfg = Config::GenerateSupercell(5, 3.4, "X", "BCC");
  supercellCfg.UpdateNeighborList(cutoffs);

  auto atomVector = cfg.GetAtomVector();
  std::set<Element> elementSet(atomVector.begin(), atomVector.end());

  PotentialEnergyEstimator peEstimator("predictor_file_WTa.json",
                                       cfg,
                                       supercellCfg,
                                       elementSet,
                                       3,
                                       3);

  elementSet.erase(Element("X"));

  VacancyMigrationBarrierPredictor barrierPredictor(cfg,
                                                      elementSet,
                                                      2,
                                                      "predictor_file_WTa.json");
 // Verificaton of dE in a loop

  auto vacancyId = cfg.GetVacancyLatticeId();
  auto initialVacancyId = vacancyId;

  size_t site1 = 15723;
  size_t site2 = 143;
  size_t site3 = 15722;

  std::cout << "Initial Vacancy ID: " << vacancyId << std::endl;

  std::pair<size_t, size_t> pair;

  // X Jump to site1
  pair = {vacancyId, site1};
  auto e1 = peEstimator.GetDe(cfg, pair);

  auto forwardBarrier1 = barrierPredictor.GetBarrier(cfg, pair);
  pair = {site1, vacancyId};  // Corrected pair for backward barrier
  auto backwardBarrier1 = barrierPredictor.GetBarrier(cfg, pair);

  std::cout << "Forward Barrier (site1) : " << forwardBarrier1 << std::endl;
  std::cout << "Backward Barrier (site1) : " << backwardBarrier1 << std::endl;

  auto e1_barrier = forwardBarrier1 - backwardBarrier1;
  std::cout << "Barrier Energy Difference (site1): " << e1_barrier << std::endl;

  cfg.LatticeJump(pair);
  std::cout << "New Vacancy ID after jump (site1): " << cfg.GetVacancyLatticeId() << std::endl;

  // X Jump to site2
  pair = {cfg.GetVacancyLatticeId(), site2};
  auto e2 = peEstimator.GetDe(cfg, pair);

  auto forwardBarrier2 = barrierPredictor.GetBarrier(cfg, pair);
  pair = {site2, cfg.GetVacancyLatticeId()};  // Corrected pair for backward barrier
  auto backwardBarrier2 = barrierPredictor.GetBarrier(cfg, pair);

  std::cout << "Forward Barrier (site2) : " << forwardBarrier2 << std::endl;
  std::cout << "Backward Barrier (site2) : " << backwardBarrier2 << std::endl;

  auto e2_barrier = forwardBarrier2 - backwardBarrier2;
  std::cout << "Barrier Energy Difference (site2): " << e2_barrier << std::endl;

  cfg.LatticeJump(pair);
  std::cout << "New Vacancy ID after jump (site2): " << cfg.GetVacancyLatticeId() << std::endl;

  // X Jump to site3
  pair = {cfg.GetVacancyLatticeId(), site3};
  auto e3 = peEstimator.GetDe(cfg, pair);

  auto forwardBarrier3 = barrierPredictor.GetBarrier(cfg, pair);
  pair = {site3, cfg.GetVacancyLatticeId()};  // Corrected pair for backward barrier
  auto backwardBarrier3 = barrierPredictor.GetBarrier(cfg, pair);

  std::cout << "Forward Barrier (site3) : " << forwardBarrier3 << std::endl;
  std::cout << "Backward Barrier (site3) : " << backwardBarrier3 << std::endl;

  auto e3_barrier = forwardBarrier3 - backwardBarrier3;
  std::cout << "Barrier Energy Difference (site3): " << e3_barrier << std::endl;

  cfg.LatticeJump(pair);
  std::cout << "New Vacancy ID after jump (site3): " << cfg.GetVacancyLatticeId() << std::endl;

  // X Jump to initial Site
  pair = {cfg.GetVacancyLatticeId(), initialVacancyId};
  auto e4 = peEstimator.GetDe(cfg, pair);

  auto forwardBarrier4 = barrierPredictor.GetBarrier(cfg, pair);
  pair = {initialVacancyId, cfg.GetVacancyLatticeId()};  // Corrected pair for backward barrier
  auto backwardBarrier4 = barrierPredictor.GetBarrier(cfg, pair);

  std::cout << "Forward Barrier (initial site) : " << forwardBarrier4 << std::endl;
  std::cout << "Backward Barrier (initial site) : " << backwardBarrier4 << std::endl;

  auto e4_barrier = forwardBarrier4 - backwardBarrier4;
  std::cout << "Barrier Energy Difference (initial site): " << e4_barrier << std::endl;

  cfg.LatticeJump(pair);
  std::cout << "New Vacancy ID after jump (initial): " << cfg.GetVacancyLatticeId() << std::endl;

  // Print and compare dE with dE from barrier
  std::cout << "Energy Change (dE) for each site jump:" << std::endl;
  std::cout << "Site 1 dE: " << e1_barrier << " , " << e1 << std::endl;
  std::cout << "Site 2 dE: " << e2_barrier << " , " << e2 << std::endl;
  std::cout << "Site 3 dE: " << e3_barrier << " , " << e3 << std::endl;
  std::cout << "Initial Site dE: " << e4_barrier << " , " << e4 << std::endl;


  std::cout << e1+e2+e3+e4 << std::endl;
  std::cout << e1_barrier + e2_barrier + e3_barrier + e4_barrier   << std::endl;





  // Verification of Barrier Prediction

  elementSet.erase(Element("X"));
  auto vacancyId = cfg.GetVacancyLatticeId();
  auto nnAtomVector = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1);

  for (auto nnAtomId : nnAtomVector)
  {


    // Jump Pair
    // { vacancyId, migratingAtomId }
    // Forward Jump

    std::cout << "Forward Jump" << std::endl;

    std::pair<size_t, size_t> JumpPair = {vacancyId,
                                          nnAtomId};

    std::cout << "{ " << JumpPair.first << "\t" << JumpPair.second << " }" << std::endl;



    std::cout << barrierPredictor.GetBarrier(cfg, JumpPair) << std::endl;


    // Backward Jump
    // { migratingAtomId, vacancyId }
    // Barrier for backward jump will be computed using the dE predicted using CE


    JumpPair = {nnAtomId,
                vacancyId};

    std::cout << "Backward Jump" << std::endl;

    std::cout << "{ " << JumpPair.first << "\t" << JumpPair.second << " }" << std::endl;


    std::cout << barrierPredictor.GetBarrier(cfg, JumpPair) << std::endl;


    std::cout << "Energy Change : " << peEstimator.GetDe(cfg,
      JumpPair) << std::endl;


    std::cout << "-----------------" << std::endl;


    break;
    }
}
*/

// #include "Home.h"
// #include "Parameter.h"

// int main(int argc, char *argv[])
// {
//   if (argc == 1) {
//     std::cout << "No input parameter filename." << std::endl;
//     return 1;
//   }
//   api::Parameter parameter(argc, argv);
//   api::Print(parameter);
//   api::Run(parameter);
// }

//
/*
void verifyDE(size_t vacId, size_t migratingAtomId, Config &cfg, VacancyMigrationPredictor vmPredictor,
PotentialEnergyEstimator pePredictor)
{
  std::cout << vacId << cfg.GetElementOfLattice(vacId)
            << "  "
            << migratingAtomId << cfg.GetElementOfLattice(migratingAtomId)
            << std::endl;

  // Forward Event
  // Atom (at migratingId) goes to Vacancy site (at VacId)

  auto barrier_De_f = vmPredictor.GetBarrierAndDiffFromLatticeIdPair(cfg,
                                                       {migratingAtomId, vacId});

  // atom -> vacancy site

  std::cout << "Forward Event" << std::endl;
  std::cout << " F Barrier : " << barrier_De_f[0];
  std::cout << " B Barrier : " << barrier_De_f[1];
  std::cout << " dE : " << barrier_De_f[2];

  std::cout << " CE dE : " << pePredictor.GetDe(cfg, {vacId, migratingAtomId});

  std::cout << std::endl;
  std::cout << cfg.GetDistanceOrder(vacId, migratingAtomId) << std::endl;


  cfg.LatticeJump({vacId, migratingAtomId});
}
*/

/*
int main() {
  // Verifying the driving force

  auto cfg = Config::ReadConfig("TiMo.POSCAR");
  cfg.UpdateNeighborList({3.20, 4.6, 5.4});

  auto supercell_cfg = Config::GenerateSupercell(4, 3.31,"Ti","BCC");
  supercell_cfg.UpdateNeighborList({3.20, 4.6, 5.4}); // for BCC Ti

  auto vacId = cfg.GetVacancyLatticeId();
  auto migratingAtomId = cfg.GetNeighborLatticeIdVectorOfLattice(vacId, 1)[0];

  VacancyMigrationPredictor vmPredictor("predictor_file.json");

  std::cout << vacId << cfg.GetElementOfLattice(vacId)
            << "  "
            << migratingAtomId << cfg.GetElementOfLattice(migratingAtomId)
            << std::endl;

  // Forward Event
  // Atom (at migratingId) goes to Vacancy site (at VacId)

  auto barrier_De_f = vmPredictor.GetBarrierAndDiffFromLatticeIdPair(cfg,
                                                       {migratingAtomId, vacId});

  auto atomVector = cfg.GetAtomVector();
  std::set<Element> elementSet(atomVector.begin(), atomVector.end());

  PotentialEnergyEstimator peEstimator("predictor_file.json",
                                       cfg,
                                       supercell_cfg,
                                       elementSet,
                                       3,
                                       3);

  // atom -> vacancy site

  std::cout << "Forward Event" << std::endl;
  std::cout << "F Barrier : " << barrier_De_f[0];
  std::cout << "B Barrier : " << barrier_De_f[1];
  std::cout << "dE : " << barrier_De_f[2];

  std::cout << std::endl;

  //cfg.LatticeJump({vacId, migratingAtomId});

  for (int i = 0; i<8; ++i)
  {
    auto tempID1 = cfg.GetNeighborLatticeIdVectorOfLattice(migratingAtomId, 1)[i];
    for (int j = 0; j<8; ++j)
    {
      auto tempID2 = cfg.GetNeighborLatticeIdVectorOfLattice(tempID1,1)[j];

      for (int z = 0; z<8; ++z)
      {
        auto tempID3 = cfg.GetNeighborLatticeIdVectorOfLattice(tempID2,1)[z];
        if (tempID3 == vacId)
        {
          if (cfg.GetElementOfLattice(tempID1) == cfg.GetElementOfLattice(tempID2) ||
              cfg.GetElementOfLattice(tempID2) == cfg.GetElementOfLattice(tempID3))
              {
                std::cout << " " << migratingAtomId << cfg.GetElementOfLattice(migratingAtomId)
                          << " " << tempID1 << cfg.GetElementOfLattice(tempID1)
                          << " " << tempID2 << cfg.GetElementOfLattice(tempID2)
                          << " " << tempID3 << cfg.GetElementOfLattice(tempID3)
                          << std::endl;
              }

        }

      }
    }


  }

  auto id_vector = {41716, 40785, 39886, 40786};

  for (auto id : id_vector)
  {
    std::cout << id << cfg.GetElementOfLattice(id) << " : ";
    for (int i = 0; i< 8; ++i)
    {
      auto temp_id = cfg.GetNeighborLatticeIdVectorOfLattice(id, 1)[i];
      std::cout << temp_id << cfg.GetElementOfLattice(temp_id) << " ";
    }
    std::cout << std::endl;
  }


  // 1
    verifyDE(vacId, migratingAtomId, cfg, vmPredictor, peEstimator);
    // 2
    verifyDE(cfg.GetVacancyLatticeId(), 39886, cfg, vmPredictor, peEstimator);

    // 3
    verifyDE(cfg.GetVacancyLatticeId(), 40786, cfg, vmPredictor, peEstimator);

    // 4
    verifyDE(cfg.GetVacancyLatticeId(), 41716, cfg, vmPredictor, peEstimator);

}
*/

/*
int main(int argc, char *argv[]) {
  // if (argc == 1) {
  //   std::cout << "No input parameter filename." << std::endl;
  //   return 1;
  // }
  // api::Parameter parameter(argc, argv);
  // api::Print(parameter);
  // api::Run(parameter);

  auto cfg = Config::ReadConfig("TiTaMoNb_Vac.POSCAR");
  cfg.UpdateNeighborList({3.20, 4.6, 5.4});


  auto old_vacancyId = cfg.GetVacancyLatticeId();

  cfg.LatticeJump({old_vacancyId, cfg.GetNeighborLatticeIdVectorOfLattice(old_vacancyId, 1)[0]});

  auto vacancyId = cfg.GetVacancyLatticeId();

  VacancyMigrationPredictor migrationPredictor("predictor_file.json");

  for (int i=0; i<8; i++){

    auto migratingAtomId = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[i];

    std::cout << "Lattice ID Pair : " << vacancyId << " " << migratingAtomId << std::endl;


    // so this GetBarrierNew function takes the lattice jump pair and will return
    // barrier for the event where the first lattice ID atom will to move to second
    // lattice ID ; one of them need to be vacancy

    // Forward Barrier
    // Migrating Atom will move to Vacancy position
    auto forward_barrier = migrationPredictor.GetBarrierNew(cfg, {migratingAtomId, vacancyId});

    auto forward_Ed = migrationPredictor.GetDiffNew(cfg, {migratingAtomId, vacancyId});

    std::cout << "Comparsion between the previous barrier and current " <<
    "barrier func for atom moving to vacany position" << std::endl;

    std::cout << "Old method barrier : " << migrationPredictor.GetBarrier(cfg, {migratingAtomId, vacancyId}) << std::endl;
    std::cout << "New method barrier : " << forward_barrier  << std::endl;

    std::cout << "Old method Ed : " << migrationPredictor.GetDiff(cfg, {migratingAtomId, vacancyId}) << std::endl;
    std::cout << "New method Ed : " << forward_Ed  << std::endl;


    // Backward Barrier
    // Assuming migrating atom and vacancy have switched places
    // So now atom is at vacancy Id
    // As our functions takes input such that specie at first Id moves to second Id
    auto backward_barrier = migrationPredictor.GetBarrierNew(cfg, {vacancyId, migratingAtomId});
    auto backward_Ed = migrationPredictor.GetDiffNew(cfg, {vacancyId, migratingAtomId});

    cfg.LatticeJump({vacancyId, migratingAtomId});

    std::cout << "Comparsion between the backward previous barrier and current " <<
    "barrier func for atom moving to vacany position" << std::endl;

    std::cout << "Old method barrier : " << migrationPredictor.GetBarrier(cfg, {migratingAtomId, vacancyId}) << std::endl;
    std::cout << "New method barrier : " << backward_barrier  << std::endl;

    std::cout << "Old method Ed : " << migrationPredictor.GetDiff(cfg, {migratingAtomId, vacancyId}) << std::endl;
    std::cout << "New method Ed : " << backward_Ed  << std::endl;

    cfg.LatticeJump({vacancyId, migratingAtomId});




    double new_forward_barrier;
    double new_forward_Ed;
    double new_backward_barrier;
    double new_backward_Ed;

    double average_barrier;
    double average_Ed;

    if (forward_barrier < backward_barrier) {

      average_Ed = (abs(forward_Ed) + abs(backward_Ed))/2;
      average_barrier = (abs(forward_barrier) + (abs(backward_barrier)-abs(backward_Ed)))/2;

      std::cout << "Average Ed : " << average_Ed << std::endl;
      std::cout << "Average barrier : " << average_barrier << std::endl;

      new_forward_barrier = average_barrier;
      new_forward_Ed = -average_Ed;
      new_backward_barrier = average_barrier + average_Ed;
      new_backward_Ed = average_Ed;
    }
    else {
      average_Ed = (abs(forward_Ed) + abs(backward_Ed))/2;
      average_barrier = ((abs(forward_barrier) - abs(forward_Ed)) + abs(backward_barrier))/2;

      std::cout << "Average Ed : " << average_Ed << std::endl;
      std::cout << "Average barrier : " << average_barrier << std::endl;
      new_forward_barrier = average_barrier + average_Ed;
      new_forward_Ed = average_Ed;
      new_backward_barrier = average_barrier ;
      new_backward_Ed = -average_Ed;

    }

    std::cout << "New barrier and Ed from the proposed solution" << std::endl;

    std::cout << "Forward Barrier : " << new_forward_barrier << std::endl;
    std::cout << "Forward Ed : " << new_forward_Ed << std::endl;
    std::cout << "Backward Barrier : " << new_backward_barrier << std::endl;
    std::cout << "Backward Ed : " << new_backward_Ed << std::endl;

    std::cout << "-----------------------" << std::endl;

    std::pair<std::pair<double, double>, std::pair<double, double>> forward_backward_info =
    {{new_forward_barrier, new_forward_Ed}, {new_backward_barrier, new_backward_Ed}};



    std::cout << "Forward Barrier // #: " << forward_backward_info.first.first << std::endl;
    std::cout << "Forward Ed // #: " << forward_backward_info.first.second << std::endl;
    std::cout << "Backward Barrier // #: " << forward_backward_info.second.first << std::endl;
    std::cout << "Backward Ed // #: " << forward_backward_info.second.second << std::endl;

    std::array<double, 3> barrier_de;
    barrier_de[0] = new_forward_barrier;
    barrier_de[1] = new_backward_barrier;
    barrier_de[2] = new_forward_Ed;

    double beta = 1 / 400 / constants::kBoltzmann;

    mc::JumpEvent event({vacancyId, migratingAtomId},
                        barrier_de,
                        beta);

    std::cout << "Defining the event " << std::endl;
    std::cout << std::endl;

    std::cout << "For Forward event " << std::endl;

    std::cout << "Forward Barrier // #: " << event.GetForwardBarrier() << std::endl;
    std::cout << "Forward Ed // #: " << event.GetEnergyChange() << std::endl;
    std::cout << "Backward Barrier: " << event.GetBackwardBarrier() << std::endl;

    auto backward_event = event.GetReverseJumpEvent();

    std::cout << "forward Barrier  for b// #: " << backward_event.GetForwardBarrier() << std::endl;
    std::cout << " Ed for backward event // #: " << backward_event.GetEnergyChange() << std::endl;
    std::cout << "backward Barrier  for b// #: " << backward_event.GetBackwardBarrier() << std::endl;

    std::cout << "**********************************************" << std::endl;

  }

}
*/

// int main(int argc, char *argv[]) {
//   // if (argc == 1) {
//   //   std::cout << "No input parameter filename." << std::endl;
//   //   return 1;
//   // }
//   // api::Parameter parameter(argc, argv);
//   // api::Print(parameter);
//   // api::Run(parameter);
//
//   auto cfg = Config::ReadConfig("TiTaMoNb_Vac.POSCAR");
//   cfg.UpdateNeighborList({3.20, 4.6, 5.4});
//
//
//   auto old_vacancyId = cfg.GetVacancyLatticeId();
//
//   cfg.LatticeJump({old_vacancyId, cfg.GetNeighborLatticeIdVectorOfLattice(old_vacancyId, 1)[0]});
//
//   auto vacancyId = cfg.GetVacancyLatticeId();
//
//
//   for (int i=0; i<8; i++){
//
//     auto migratingAtomId = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[i];
//
//     std::cout << "Lattice ID Pair : " << vacancyId << " " << migratingAtomId << std::endl;
//
//     VacancyMigrationPredictor migrationPredictor("predictor_file.json");
//
//     auto forward_barrier = migrationPredictor.GetBarrier(cfg, {vacancyId, migratingAtomId});
//     auto forward_Ed = migrationPredictor.GetDiff(cfg, {vacancyId, migratingAtomId});
//
//     std::cout << migratingAtomId << cfg.GetElementOfLattice(migratingAtomId) << " ";
//     std::cout << forward_barrier << " " << forward_Ed << std::endl;
//
//     cfg.LatticeJump({vacancyId, migratingAtomId});
//
//     // how to get the barrier and Ed without swapping the positions
//
//
//     auto backward_barrier = migrationPredictor.GetBarrier(cfg, {vacancyId, migratingAtomId});
//     auto backward_Ed = migrationPredictor.GetDiff(cfg, {vacancyId, migratingAtomId});
//
//     std::cout << migratingAtomId << cfg.GetElementOfLattice(migratingAtomId) << " ";
//     std::cout << backward_barrier << " " << backward_Ed << std::endl;
//
//     cfg.LatticeJump({vacancyId, migratingAtomId});
//
//
//
//
//     double new_forward_barrier;
//     double new_forward_Ed;
//     double new_backward_barrier;
//     double new_backward_Ed;
//
//     double average_barrier;
//     double average_Ed;
//
//     if (forward_barrier < backward_barrier) {
//
//       average_Ed = (abs(forward_Ed) + abs(backward_Ed))/2;
//       average_barrier = (abs(forward_barrier) + (abs(backward_barrier)-abs(backward_Ed)))/2;
//
//       std::cout << "Average Ed : " << average_Ed << std::endl;
//       std::cout << "Average barrier : " << average_barrier << std::endl;
//       new_forward_barrier = average_barrier;
//       new_forward_Ed = -average_Ed;
//       new_backward_barrier = average_barrier + average_Ed;
//       new_backward_Ed = average_Ed;
//     }
//     else {
//       average_Ed = (abs(forward_Ed) + abs(backward_Ed))/2;
//       average_barrier = ((abs(forward_barrier) - abs(forward_Ed)) + abs(backward_barrier))/2;
//
//       std::cout << "Average Ed : " << average_Ed << std::endl;
//       std::cout << "Average barrier : " << average_barrier << std::endl;
//       new_forward_barrier = average_barrier + average_Ed;
//       new_forward_Ed = average_Ed;
//       new_backward_barrier = average_barrier ;
//       new_backward_Ed = -average_Ed;
//
//     }
//
//     std::cout << "New barrier and Ed from the proposed solution" << std::endl;
//
//     std::cout << "Forward Barrier : " << new_forward_barrier << std::endl;
//     std::cout << "Forward Ed : " << new_forward_Ed << std::endl;
//     std::cout << "Backward Barrier : " << new_backward_barrier << std::endl;
//     std::cout << "Backward Ed : " << new_backward_Ed << std::endl;
//
//     std::cout << "-----------------------" << std::endl;
//
//     std::pair<std::pair<double, double>, std::pair<double, double>> forward_backward_info =
//     {{new_forward_barrier, new_forward_Ed}, {new_backward_barrier, new_backward_Ed}};
//
//
//
//     std::cout << "Forward Barrier // #: " << forward_backward_info.first.first << std::endl;
//     std::cout << "Forward Ed // #: " << forward_backward_info.first.second << std::endl;
//     std::cout << "Backward Barrier // #: " << forward_backward_info.second.first << std::endl;
//     std::cout << "Backward Ed // #: " << forward_backward_info.second.second << std::endl;
//
//     std::array<double, 3> barrier_de;
//     barrier_de[0] = new_forward_barrier;
//     barrier_de[1] = new_backward_barrier;
//     barrier_de[2] = new_forward_Ed;
//
//     double beta = 1 / 400 / constants::kBoltzmann;
//
//     mc::JumpEvent event({vacancyId, migratingAtomId},
//                         barrier_de,
//                         beta);
//
//     std::cout << "Defining the event " << std::endl;
//     std::cout << std::endl;
//
//     std::cout << "For Forward event " << std::endl;
//
//     std::cout << "Forward Barrier // #: " << event.GetForwardBarrier() << std::endl;
//     std::cout << "Forward Ed // #: " << event.GetEnergyChange() << std::endl;
//     std::cout << "Backward Barrier: " << event.GetBackwardBarrier() << std::endl;
//
//     auto backward_event = event.GetReverseJumpEvent();
//
//     std::cout << "forward Barrier  for b// #: " << backward_event.GetForwardBarrier() << std::endl;
//     std::cout << " Ed for backward event // #: " << backward_event.GetEnergyChange() << std::endl;
//     std::cout << "backward Barrier  for b// #: " << backward_event.GetBackwardBarrier() << std::endl;
//
//   }
//
// }
//
//// int main() {
////   mc::ThermodynamicAveraging thermodynamic_avg(0);
////   thermodynamic_avg.AddEnergy(3838);
////   thermodynamic_avg.AddEnergy(3838);
////   thermodynamic_avg.AddEnergy(3838);
//   thermodynamic_avg.AddEnergy(3838);
//
//   thermodynamic_avg.GetThermodynamicAverage(1);
// }

// double GetGeometricMeanMain(std::vector<double> xSv_vector, double z = 1.0) {
//   double product = 1;
//   double n = static_cast<double>(xSv_vector.size());
//
//   for (auto& xSv : xSv_vector) {
//     product *= xSv;
//   }
//
//   return pow(product, z/n);
// }
//
// double GetxSvMain(const Element migrating_ele) {
//   return migrating_ele.GetElectronegativity() * migrating_ele.GetSv();
// }
//
// int main(int argc, char *argv[]){
//
//   auto start = std::chrono::high_resolution_clock::now();
//
// Call the function to measure
// Config cfg = Config::ReadPoscar("Ti_Ni_BCC_Vacancy.poscar");
// cfg.UpdateNeighborList({3.20,4.6,5.4}); // for BCC Ti_Ni
//
// auto supercell_cfg = Config::GenerateSupercell(4, 3.31,"Ti","BCC");
// supercell_cfg.UpdateNeighborList({3.20, 4.6, 5.4}); // for BCC Ti

//////////////////////// FCC System ///////////////////

//     auto cfg_Al = Config::ReadCfg("start.cfg");
//     cfg_Al.UpdateNeighborList({3.5, 4.8, 5.3});
//     std::cout << "finish reading config" << std::endl;
//     //auto temp = FindAllLatticeClustersMain(cfg,3,3,{});
//
//
//
//      //auto cfg = Config::ReadPoscar("Al.poscar");
//      auto supercell_cfg = Config::GenerateSupercell(4,4.046,"Al","FCC");
//      supercell_cfg.UpdateNeighborList({3.5, 4.8, 5.3});
//
//
//     auto atomVector = cfg_Al.GetAtomVector();
//     std::set<Element> element_set(atomVector.begin(), atomVector.end());
//     for (const auto& element : element_set) {
//           std::cout << element.GetElementString() << " ";
//       }
//
//     std::cout << std::endl;
//
//     // first nearest neighbours list
//     auto nn_list1 = cfg_Al.GetNeighborLatticeIdVectorOfLattice(1556,1);
//
//     std::cout << "1st Neighbours List: { ";
//     for(auto id : nn_list1){
//       std::cout << id << ", ";
//     }
//     std::cout << "}" << std::endl;
//
//   std::string predictor_file_name = "predictor_file.json";
//
//
//     PotentialEnergyEstimator pe_estimator(predictor_file_name,
//                                           cfg_Al,
//                                           supercell_cfg,
//                                           element_set,
//                                           3,
//                                           3);
//
//    auto de = pe_estimator.GetDe(cfg_Al, {0, 647});
//    std::cout << "0 -> 647 : " << de << std::endl;
//    std::cout << "647 -> 0 : " << pe_estimator.GetDe(cfg_Al, {647, 0}) << std::endl;
//
//
//     // first nearest neighbours list
//     auto nn_list2 = cfg.GetNeighborLatticeIdVectorOfLattice(1556,2);
//
//     std::cout << "2nd Neighbours List: { ";
//     for(auto id : nn_list2){
//       std::cout << id << ", ";
//     }
//     std::cout << "}" << std::endl;
//
//     // first nearest neighbours list
//     auto nn_list3 = cfg.GetNeighborLatticeIdVectorOfLattice(1556,3);
//
//     std::cout << "3rd Neighbours List: { ";
//     for(auto id : nn_list3){
//       std::cout << id << ", ";
//     }
//     std::cout << "}" << std::endl;
//
//
//
//     std::pair<size_t,size_t> lattice_id_pair = {19072,19633};
//
//     size_t max_bond_order = 3;
//
//     auto nn = cfg.GetNeighboringLatticeIdSetOfPair(lattice_id_pair, max_bond_order);
//
//     std::cout << "Nearest Neighbours of the lattice ID Pair: { ";
//     for (const auto id : nn){
//       std::cout << id << ", ";
//     }
//     std::cout << "}" << std::endl;
//
//     std::cout << "Size of the nn set : " << nn.size() << std::endl;
//
//     auto center = cfg.GetLatticePairCenter(lattice_id_pair);
//
//     std::cout << "Center of the lattice Id pair : " << center << std::endl;
//
//     auto matrix = cfg.GetLatticePairRotationMatrix( lattice_id_pair);
//
//     std::cout << "Lattice Pair Rotation Matrix: " << matrix << std::endl;
//
//     auto ss = GetSymmetricallySortedLatticeVectorMMM(cfg, lattice_id_pair, max_bond_order);
//
//     std::cout << "SymmetricallySortedLatticeVectorMMM: {";
//     for(auto id : ss) {
//       std::cout << id << " , " << cfg.GetRelativePositionOfLattice(id).transpose() << std::endl;
//     }
//     std::cout << "Size of ss: " << ss.size() << std::endl;
//
//     auto ss_mm2 = GetSymmetricallySortedLatticeVectorMM2(cfg, lattice_id_pair, max_bond_order);
//
//     std::cout << "SymmetricallySortedLatticeVectorMM2: {";
//     for(auto id : ss_mm2) {
//       std::cout << id << " , " << cfg.GetRelativePositionOfLattice(id).transpose() << std::endl;
//
//     }
//
//
//     /// For getting average parameter mapping from lattice cluster MMM
//
//     // Lattice Id Jump Pair
//     std::pair<size_t, size_t>
//     lattice_id_jump_pair = {0, cfg.GetNeighborLists()[0][0][0]};
//
//     std::cout << "Lattice Id Jump Pair MMM : {" << lattice_id_jump_pair.first
//       << ", " << lattice_id_jump_pair.second << "}" << std::endl;
//
//     auto lattice_vector = GetSymmetricallySortedLatticeVectorMMM(cfg, lattice_id_jump_pair, max_bond_order);
//
//     std::vector<std::pair<size_t, size_t>> lattice_vector_hashmap;
//     // original lattice id, index corresponding to MMM symmetry
//
//
//
//     std::vector<std::pair<size_t, Eigen::RowVector3d>> singlet_vector;
//     std::unordered_set<std::vector<size_t>, boost::hash<std::vector<size_t>>> first_pair_set;
//     std::unordered_set<std::vector<size_t>, boost::hash<std::vector<size_t>>> second_pair_set;
//     std::unordered_set<std::vector<size_t>, boost::hash<std::vector<size_t>>> third_pair_set;
//
//
//     for (size_t id1 = 0; id1 < lattice_vector.size(); ++id1){
//       Eigen::RowVector3d pos = cfg.GetRelativePositionOfLattice(id1).transpose();
//       singlet_vector.push_back({id1, pos});
//       for (size_t id2 = 0; id2<id1; id2++){
//         switch (cfg.GetDistanceOrder(lattice_vector[id1], lattice_vector[id2]))
//         {
//         case 1: first_pair_set.emplace(std::vector<size_t>{id1, id2});
//           break;
//
//         case 2: second_pair_set.emplace(std::vector<size_t>{id1, id2});
//           break;
//
//         case 3: third_pair_set.emplace(std::vector<size_t>{id1, id2});
//           break;
//
//         }
//       }
//     }
//
//
//     for(auto cluster : singlet_vector){
//       std::cout << cluster.first << " : " << cluster.second << std::endl;
//     }
//
//     // std::sort(singlet_vector.begin(), singlet_vector.end(), IsClusterSmallerSymmetricallyMMM);
//
//     for(auto cluster : singlet_vector){
//       std::cout << cluster.first << " : " << cluster.second << std::endl;
//     }
//
//

// PotentialEnergyEstimator
// estimator("quartic_coefficients.json",cfg,supercell_cfg,
//            element_set,3,3);
//
// auto total_energy = estimator.GetEnergy(cfg);

// mc::CanonicalMcSerial mc_simulation(cfg,supercell_cfg,100,100,10,10,0,0,500,element_set,3,3,"quartic_coefficients.json");

//    auto lattice_id_pair = {13613,805};
//
//    auto nn_list = cfg.GetNeighborLatticeIdVectorOfLattice(13613,1);

//    auto temp = getUniqueNeighbors(cfg, 19072,19633);
//
//    std::cout << "nn upto 3rd nn : { ";
//    for(const auto& id : temp){
//    std::cout << id << ", ";
//    }
//    std::cout << "}" << std::endl;
//
//    std::cout << temp.size() << std::endl;

//    for(auto& id : temp){
//      FindPointGroupSymmetry(cfg,{id});
//    }

//    std::vector<size_t> temp_vector(temp.begin(), temp.end());
//
//    FindPointGroupSymmetry(cfg, temp_vector);

// auto sorted_ids = SymmetricallySortLatticeIDs(cfg,temp_vector);

// auto temp = FindAllLatticeClusters(cfg,2,3,{1382,805});

// for(auto& cluster : temp){
//   std::vector<size_t> lattice_cluster = cluster.GetLatticeIdVector();
//   FindPointGroupSymmetry(cfg,lattice_cluster);
// }

// Convert the unordered_set to a vector for indexed access
//    std::vector<size_t> temp_vector(temp.begin(), temp.end());
//
//    // Generate all unique pairs and pass them to FindPointGroupSymmetry
//    std::vector<std::vector<size_t>> pair_cluster_vector;
//    for (size_t i = 0; i < temp_vector.size(); i++) {
//        for (size_t j = i + 1; j < temp_vector.size(); j++) {
//            // Create a pair cluster and pass it to FindPointGroupSymmetry
//            std::vector<size_t> cluster = {temp_vector[i], temp_vector[j]};
//            pair_cluster_vector.push_back(cluster);
//        }
//    }

// FindEquivalentClusters(cfg,pair_cluster_vector);
//
//
//   for(auto& id1 : nn_list){
//    std::cout << id1 << cfg.GetElementOfLattice(id1) << ", ";
//   }
//

//  std::pair<size_t,size_t> lattice_id_pair = {19072,19633};

// std::cout << "Vacancy Lattice ID : " << cfg.GetVacancyLatticeId() << std::endl;
//
// auto center = GetLatticePairCenter(cfg,lattice_id_pair);
//
// std::cout << center << std::endl;

//  Eigen::Vector3d rdv = cfg.GetRelativeDistanceVectorLattice(19072,19633);
//
//  std::cout << rdv << std::endl;
//
//  Eigen::Matrix3d basis = cfg.GetBasis();
//
//  Eigen::Vector3d result = basis * rdv; // Matrix-vector multiplication
//  auto result1 = result.normalized();
//
//  std::cout << "Result of multiplication: " << result1 << std::endl;

// std::cout << cfg.GetRelativeDistanceVectorLattice(lattice_id_pair.first, lattice_id_pair.second) << std::endl;
//
// std::cout << cfg.GetBasis() << std::endl;
//
// std::cout << (cfg.GetRelativeDistanceVectorLattice(lattice_id_pair.first, lattice_id_pair.second).transpose() * cfg.GetBasis()).normalized() << std::endl;
// std::cout << (cfg.GetRelativeDistanceVectorLattice(lattice_id_pair.first, lattice_id_pair.second).transpose() * cfg.GetBasis()).normalize() << std::endl;

// auto rot_mat = GetLatticePairRotationMatrix(cfg,lattice_id_pair);
//
// std::cout << rot_mat << std::endl;

//  VacancyMigrationPredictor barrier_predictor("quartic_coefficients.json",cfg, element_set);
//
//  // auto barrier_de = barrier_predictor.GetBarrierAndDiffFromLatticeIdPair(cfg,lattice_id_pair);
//
//  // std::cout << "Barrier: " << barrier_de.first << "; De: " << barrier_de.second << std::endl;
//
//  std::string json_filename = "quartic_coefficients.json";
//  std::string temp_filename = "";
//  Eigen::RowVector3d vacancy_trajectory = {0,0,0};
//  std::pair<double, double> barrier_and_dE = {1.2, 0.12};
//
//  double beta = 1/(constants::kBoltzmann * 600);

//  mc::JumpEvent jump_event(lattice_id_pair, barrier_and_dE, beta);

//  std::cout << "Forward Rate: " << jump_event.GetForwardRate() << std::endl;

// mc::KineticMcFirstAbstract kmc_first(cfg, supercell_cfg, 10, 10, 10, 1, 0, 0, 0, 400, element_set, 3, 3, json_filename, temp_filename, false, vacancy_trajectory);

// mc::KineticMcChainOmpi kmc_ompi(cfg, supercell_cfg, 10, 10, 10, 1, 0, 0, 0, 400, element_set, 3, 3, json_filename, temp_filename, false, vacancy_trajectory);

//  auto temp = cfg.GetRelativeDistanceVectorLattice(980,981);
//
//  // PotentialEnergyEstimator estimator(json_filename,cfg, supercell_cfg, element_set, 3, 3);
//
//  auto clsuters = FindAllLatticeClusters(supercell_cfg, 3,3, {});

// mc::KineticMcFirstMpi kmc_mpi(cfg, supercell_cfg, 10,
//                                  10,10,1,0,0,0,1000,element_set,
//                                  3,3,json_filename, temp_filename,
//                                  false, vacancy_trajectory);

//  PotentialEnergyEstimator estimator(json_filename, cfg, supercell_cfg, element_set, 3, 3);
//
//  auto cfg_new = Config::ReadConfig("Ti_Ni_BCC_Vacancy.poscar");
//  cfg_new.UpdateNeighborList({3.5, 4.8, 5.3});

// auto kmc = mc::KineticMcChainOmpi(cfg,10,10,10,1,0,0,0,400,element_set, json_filename, temp_filename, false, vacancy_trajectory);

// auto end1 = std::chrono::high_resolution_clock::now();
//
// auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start).count();
//
// std::cout << "Time taken for mc_simulation to initialize: " << duration1 << " milliseconds" << std::endl;

// mc_simulation.Simulate();

// auto temp = IdentifyLatticeClusterType(cfg,{1379,1650,1072});
// std::cout << temp << std::endl;
// cfg.GetDistanceOrder(1379,1650);
// cfg.GetDistanceOrder(1072,1650);
// cfg.GetDistanceOrder(1379,1072);

// printNestedList(n_list_lattice);

// std::cout << "Distance Order: " << cfg.GetDistanceOrder(1090,1667) << std::endl;

// PotentialEnergyEstimator
// estimator("quartic_coefficients.json",cfg,supercell_cfg,
//           element_set,3,3);
//
// auto temp = estimator.GetEnergy(cfg);
// std::cout << temp << std::endl;
// std::cout << cfg.GetNumLattices() << std::endl;

// End timing

///////////////////////// BCC System ////////////////////////

//  Config start_cfg = Config::ReadConfig("TiTaMoNb_Vac.POSCAR");
//  start_cfg.UpdateNeighborList({3.20,4.6,5.4}); // for Ti as the base matrix
////
//  auto atomVector = start_cfg.GetAtomVector();
//  std::set<Element> element_set(atomVector.begin(), atomVector.end());
//  ShortRangeOrder sro_start(start_cfg, element_set);
//
//  size_t shell_number = 1;
//
//  auto start_wcp = sro_start.FindWarrenCowley(shell_number);
//
//  Config end_cfg = Config::ReadConfig("end_NbTaTiW_800.cfg");
//  end_cfg.UpdateNeighborList({3.20,4.6,5.4}); // for Ti as the base matrix
//
//  ShortRangeOrder sro_end(end_cfg, element_set);
//  auto end_wcp = sro_end.FindWarrenCowley(shell_number);
//
//  std::cout << std::endl;
//  for (auto ele_pair : start_wcp) {
//    std::cout << ele_pair.first << " : " << ele_pair.second << std::endl;
//  }
//
//  std::cout << std::endl;
//  for (auto ele_pair : end_wcp) {
//    std::cout << ele_pair.first << " : " << ele_pair.second << std::endl;
//  }
//

// sro param without considering vacancy in the element set
//  ShortRangeOrder sro(start_cfg, element_set);
//  auto wcp_map = sro.FindWarrenCowley(1);
//  for (auto wcp : wcp_map) {
//    std::cout << wcp.first << " : " << wcp.second << std::endl;
//  }

//  api::Parameter parameter(argc, argv);
//  api::Print(parameter);
//  api::Run(parameter);

///////////////////// Testing Second Order KMC Algo ///////////////////////////
//  VacancyMigrationPredictor migrationPredictor("predictor_file.json");
//
//  auto vacancyId = start_cfg.GetVacancyLatticeId();
//  std::cout << "Vacancy Lattice Id : " << vacancyId << std::endl;
//
//
//  auto firstNN = start_cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1);
//
//  size_t previous_lattice_id{};
//
//  for (auto neighbouringId : firstNN) {
//
//    if (previous_lattice_id)
//    {
//    // putting back vacancy to its original place
//    start_cfg.LatticeJump({vacancyId, previous_lattice_id});
//    }
//    std::cout << "Initial vacancy position : " << vacancyId << std::endl;
//    auto neighbourIDVector = start_cfg.GetNeighborLatticeIdVectorOfLattice(neighbouringId, 1);
//    std::cout << "Id1_vac : Id2 : forwardBarrier : dE : backwardBarrier : -dE" <<
//    ": backwardBarrier_xSv : backward_dE : id2_vac" << std::endl;
//
//    // move vacancy to neighbouringId
//    start_cfg.LatticeJump({vacancyId, neighbouringId});
//    previous_lattice_id = neighbouringId;
//
//    for (auto id : neighbourIDVector) {
//      std::cout << neighbouringId << start_cfg.GetElementOfLattice(neighbouringId)
//      << " : " << id << start_cfg.GetElementOfLattice(id) << " : ";
//      auto forwardBarrier = migrationPredictor.GetBarrier(start_cfg, {neighbouringId, id});
//      auto dE = migrationPredictor.GetDiff(start_cfg, {neighbouringId, id});
//      auto backwardBarrier = forwardBarrier - dE;
//
//      std::cout << forwardBarrier << " : " << dE << " : " << backwardBarrier <<
//      " : " << -dE << " : ";
//
//      // Backward barrier using xSv
//      // neighbouringId has the vacancy
//
//      start_cfg.LatticeJump({neighbouringId,id});
//
//      auto backwardBarrier_xSv = migrationPredictor.GetBarrier(start_cfg, {neighbouringId,id});
//      auto backward_dE = migrationPredictor.GetDiff(start_cfg, {neighbouringId,id});
//
//      std::cout << backwardBarrier_xSv << " : " << backward_dE << " : " << id <<
//      start_cfg.GetElementOfLattice(id) << std::endl;
//      // std::cout << "Verificatin : " << id << start_cfg.GetElementOfLattice(id) << std::endl;
//
//      start_cfg.LatticeJump({neighbouringId,id});
//      // std::cout << "Verificatin : " << id << start_cfg.GetElementOfLattice(id) << std::endl;
//
//
//    }
//    std::cout << "----------------------------------------------------" << std::endl;
//  }

////////////////////////// Verification of GetDiff function ////////////////////

//  auto vac_id = start_cfg.GetVacancyLatticeId();
//  auto neighbouring_id = start_cfg.GetNeighborLatticeIdVectorOfLattice(vac_id, 1)[0];

//  start_cfg.LatticeJump({vac_id, neighbouring_id});
//  VacancyMigrationPredictor migrationPredictor("predictor_file.json");
//
//  auto dE = migrationPredictor.GetDiff(start_cfg, {vac_id, neighbouring_id});
//
//  std::cout << "dE : " << dE << std::endl;
//
//  auto forward_barrier = migrationPredictor.GetBarrier(start_cfg, {vac_id, neighbouring_id});
//
//  start_cfg.LatticeJump({vac_id, neighbouring_id});
//
//  auto backward_barrier = migrationPredictor.GetBarrier(start_cfg, {vac_id, neighbouring_id});
//
//
//  std::cout << "Forward Barrier : " << forward_barrier << std::endl;
//  std::cout << "Backward Barrier : " << backward_barrier << std::endl;
//
//  std::cout << migrationPredictor.GetDiff(start_cfg, {vac_id, neighbouring_id}) << std::ednl;
//  start_cfg.LatticeJump({vac_id, neighbouring_id});

//
//  auto vac_id = cfg.GetVacancyLatticeId();
//
//  std::cout << "Vacancy Lattice ID: " << vac_id << std::endl;
//
//  auto vac_neighbours = cfg.GetNeighborLatticeIdVectorOfLattice(vac_id, 1);
//
//  std::string predictor_file = "predictro";
//  VacancyMigrationPredictor mig_predictor(predictor_file);
//
//  std::cout << "Neighbours : " << std::endl;
//  for (auto id : vac_neighbours) {
//    std::cout << id << " ; " << mig_predictor.GetDiff(cfg, {vac_id, id}) <<
//    " ; " << mig_predictor.GetBarrier(cfg, {vac_id, id}) << std::endl;
//  }
//
//  std::pair<size_t, size_t> lattice_id_pair = {6939, vac_id};
//
//  auto formFunc = mig_predictor.GetBarrierAndDiffFromLatticeIdPair(cfg, lattice_id_pair);
//  std::cout << formFunc.first << "; " << formFunc.second << std::endl;
//
//
//  std::cout << mig_predictor.GetBarrier(cfg, lattice_id_pair) << "; " <<
//  mig_predictor.GetDiff(cfg, lattice_id_pair) << std::endl;
//
//  std::cout << "{ " ;
//  for (auto ele : cfg.GetAtomVector()) {
//    std::cout << ele.GetElementString() << " ,";
//  }
//  std::cout << std::endl;
//
//  std::unordered_map<std::string, int> counts;
//
//    // Count occurrences
//    for (const auto& item : cfg.GetAtomVector()) {
//        counts[item.GetElementString()]++;
//    }
//
//    // Print the counts
//    std::cout << "Occurrences:" << std::endl;
//    for (const auto& pair : counts) {
//        std::cout << pair.first << ": " << pair.second << std::endl;
//    }
//
// auto conc = cfg.GetConcentration();
// for (auto ele : conc) {
//  std::cout << ele.first.GetElementString() << " : " << ele.second << std::endl;
// }
//
//  auto nn_list_1 = cfg.GetNeighborLatticeIdVectorOfLattice(0,1);
//  auto nn_list_2 = cfg.GetNeighborLatticeIdVectorOfLattice(0,2);
//  auto nn_list_3 = cfg.GetNeighborLatticeIdVectorOfLattice(0,3);
//
//  std::cout << "First NN : " << nn_list_1.size() << std::endl;
//  std::cout << "Second NN : " << nn_list_2.size() - nn_list_1.size() << std::endl;
//  std::cout << "Third NN : " << nn_list_3.size()- nn_list_2.size() << std::endl;
//
//
//  std::cout << "Lattice ID 100 : " << cfg.GetAtomVector()[100].GetElementString()
//   << " : " << cfg.GetElementOfAtom(100) << std::endl;
//
//   std::cout << "Lattice ID 1000 : " << cfg.GetAtomVector()[1000].GetElementString()
//   << " : " << cfg.GetElementOfAtom(1000) << std::endl;
//
//   std::cout << "Lattice ID 500 : " << cfg.GetAtomVector()[500].GetElementString()
//   << " : " << cfg.GetElementOfAtom(500) << std::endl;
//
//   std::cout << "atom neighour of atom id 100 : " << std::endl;
//   for (auto id : cfg.GetNeighborAtomIdVectorOfAtom(100, 1)){
//    std::cout << id << " : " << cfg.GetElementOfAtom(id) << std::endl;
//   }
//
//   std::cout << "lattice neighour of lattice id 100 : " << std::endl;
//   for (auto id : cfg.GetNeighborLatticeIdVectorOfLattice(100, 1)){
//    std::cout << id << " : " << cfg.GetElementOfLattice(id) <<  std::endl;
//   }
//
//  auto atomVector = cfg.GetAtomVector();
//  std::set<Element> element_set(atomVector.begin(), atomVector.end());
//
//  for (const auto& element : element_set) {
//    std::cout << element.GetElementString() << " ";
//  }
//  std::cout << std::endl;
//
//  ShortRangeOrder sro(cfg, element_set);
//  auto wcp = sro.FindWarrenCowley(1);
////
//  for (auto wc : wcp ) {
//    std::cout << wc.first << " : " << wc.second << std::endl;
//  }

// Verification of warren cowley parameter

//  auto cfg_Al = Config::ReadCfg("start.cfg");
//  cfg_Al.UpdateNeighborList({3.5, 4.8, 5.3});
//
//  auto atom_vector_start = cfg_Al.GetAtomVector();
//  std::set<Element> element_set_start(atom_vector_start.begin(), atom_vector_start.end());
//
//  std::cout << "Element Set for start.cfg : ";
//  for (auto ele : element_set_start) {
//    std::cout << ele.GetElementString() << " ,";
//  }
//  std::cout << std::endl;
//
//  ShortRangeOrder sro_start(cfg_Al, element_set_start);
//
//  std::cout << "sro_start has been declared" << std::endl;
//  auto wcp_start = sro_start.FindWarrenCowley(1);
//
//  for (auto wcp : wcp_start) {
//    std::cout << wcp.first << " : " << wcp.second << std::endl;
//  }

//  std::string predictor_file_name = "predictor_file.json";
//
//
//  VacancyMigrationPredictor vac_predictor(predictor_file_name);

// Second Order KMC code

// }
//
//  Eigen::RowVector3d cartesian_pos = cfg.GetCartesianPositionOfLattice(1000);
//  Eigen::RowVector3d cartesian_pos711 = cfg.GetCartesianPositionOfLattice(711);
//
// auto nn_711 = cfg.GetNeighborLatticeIdVectorOfLattice(711, 1);
//
//  std::cout << "NN of 711 : " << std::endl;
//
//  for (auto& id : nn_711){
//    std::cout << id << std::endl;
//  }
//
//  std::cout << "Cartesian Position of 1000 : " << cartesian_pos << std::endl;
//  std::cout << "Cartesian Position of 711 : " << cartesian_pos711 << std::endl;
//
//
//
//
//    std::cout << "********************** xSv method **********************" << std::endl;
//
//    std::cout << " Lattice ID Jump Pair : { " << lattice_id_jump_pair.first <<
//                ", " << lattice_id_jump_pair.second << " }" << std::endl;
//
//    // Combined 1st NN of lattice Id Jump Pair
//    auto nn_1 = cfg.GetNeighboringLatticeIdSetOfPair(lattice_id_jump_pair, 1); // first nearest neighbours
//
//    std::cout << "First Nearest Neighbours : " << std::endl;
//
//    for (auto& id : nn_1) {
//      Eigen::RowVector3d cartesian_pos = cfg.GetCartesianPositionOfLattice(id).transpose();
//      std::cout << id << " :  { " << cartesian_pos << " }" << std::endl;
//    }
//
//
//  std::unordered_map<std::string, size_t> Sv_map = {
//    {"W", 6},
//    {"Mo", 6},
//    {"Ta", 5},
//    {"Nb", 5},
//    {"V", 5},
//    {"Hf", 4},
//    {"Zr", 4},
//    {"Ti", 4},
//    {"Al", 3},
//    {"Ni", 10}
//  };
//
//
//  // Assuming
//  // Initial Position = 1000
//  // Final Position = 711
//  std::vector<size_t> i1_positions; // I1
//  std::vector<size_t> i2_positions; // I2
//  std::vector<size_t> f1_positions; // F1
//  std::vector<size_t> f2_positions; // F2
//
//
//  for (auto& id : nn_1) {
//
//
//
//    Eigen::RowVector3d cartesian_pos_id = cfg.GetCartesianPositionOfLattice(id).transpose();
//
//    // std::cout << id << " ; { " << cartesian_pos_id << " } ; ";
//
//    auto ele = cfg.GetElementOfLattice(id);
//
//    std::cout << id << " : " << ele.GetElementString() << " : " << Sv_map.at(ele.GetElementString()) * ele.GetElectronegativity()<< " : ";
//
//    auto initial_BO = cfg.GetDistanceOrder(lattice_id_jump_pair.first, id);
//    auto final_BO = cfg.GetDistanceOrder(lattice_id_jump_pair.second, id);
//
//    // condition for I1 atoms
//    if (initial_BO == 1 && final_BO == 2) {
//      i1_positions.push_back(id);
//      std::cout << "I1 atom" << std::endl;
//    }
//
//    // condition for F1 atoms
//    else if (initial_BO == 2 && final_BO == 1) {
//      f1_positions.push_back(id);
//      std::cout << "F1 atom" << std::endl;
//    }
//
//    // condition for I2 atoms
//    else if (initial_BO == 1 && (final_BO >= 3)) {
//      i2_positions.push_back(id);
//      std::cout << "I2 atom" << std::endl;
//    }
//
//    // condition for F2 atoms
//    else if ((initial_BO >= 3) && final_BO == 1 ) {
//      f2_positions.push_back(id);
//      std::cout << "F2 atom" << std::endl;
//    }
//
//  }
//
//  // verification of electronegativity values
//
//  std::vector<std::string> ele_vector = {"Al","Hf","Zr","Ti","Ta","Nb","V","Mo","W"};
//
//
//
//
//  for (auto& ele_string : ele_vector) {
//    Element ele(ele_string);
//    std::cout << ele.GetElementString() << " : " << Sv_map.at(ele_string) * ele.GetElectronegativity() << std::endl;
//  }
//
//
//  std::vector<double> temp_ = {1, 2, 3, 4};
//  auto result = GetGeometricMeanMain(temp_, 4);
//
//  std::cout << result << std::endl;
//
//  auto ele_check = cfg.GetElementOfLattice(1307);
//
//  std::cout << 1307 << " : " << ele_check.GetElementString() << " : " << Sv_map.at(ele_check.GetElementString()) * ele_check.GetElectronegativity()<< std::endl;
//
//
//  std::cout << ele_check.GetSv() << std::endl;
//
//
//  std::cout << GetxSvMain(ele_check) << std::endl;
//
//  std::string file_name = "predictor";
//
//  std::set<Element> element_null_set{};
//
//
//  VacancyMigrationPredictor eb_predictor(file_name);
//
//  eb_predictor.GetBarrier(cfg, lattice_id_jump_pair);
//
//
//  std::cout << GetDf(cfg, lattice_id_jump_pair.first) << std::endl;
//  std::cout << GetDf(cfg, lattice_id_jump_pair.second) << std::endl;
//
//  std::cout << GetDf(cfg, lattice_id_jump_pair.first) - GetDf(cfg, lattice_id_jump_pair.second) << std::endl;
//
//
//  std::cout << eb_predictor.GetDiff(cfg, lattice_id_jump_pair);
//
//  std::string predictor_file_name = "predictor_file.json";
//  auto atomVector = cfg.GetAtomVector();
//  std::set<Element> element_set(atomVector.begin(), atomVector.end());
//  for (const auto& element : element_set) {
//    std::cout << element.GetElementString() << " ";
//  }
//
//  std::cout << lattice_id_jump_pair.first << cfg.GetElementOfLattice(lattice_id_jump_pair.first) <<
//  " : " << lattice_id_jump_pair.second << cfg.GetElementOfLattice(lattice_id_jump_pair.second) << std::endl;
//
//
//  PotentialEnergyEstimator pe_estimator(predictor_file_name, cfg, supercell_cfg, element_set, 3, 3);
//
//  std::cout << pe_estimator.GetDe(cfg, lattice_id_jump_pair) << std::endl;
//
//  std::cout << lattice_id_jump_pair.first << cfg.GetElementOfLattice(lattice_id_jump_pair.first) <<
//  " : " << lattice_id_jump_pair.second << cfg.GetElementOfLattice(lattice_id_jump_pair.second) << std::endl;
//
//
//
//
//
//
//
//
//
//
//

//
//
//    auto end = std::chrono::high_resolution_clock::now();
//
//    // Calculate elapsed time in milliseconds
//    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
//
//    size_t max_bond_order = 3;
//    cfg.WriteLattice("lattice_bcc_Ti_Ni.txt",max_bond_order);
//
//    std::cout << "Time taken for Simulation: " << duration << " milliseconds" << std::endl;

// }

//
//
//
//  Config config;
//  if (parameter.map_filename_.empty()) {
//    // Need to make a generalized function
//    // config = Config::ReadPoscar(parameter.config_filename_);
//    config = Config::ReadCfg(parameter.config_filename_);
//
//    config.UpdateNeighborList(parameter.cutoffs_);
//    //std::cout << "Read the Configuration" << std::endl;
//
//  } //else {
//    //config = Config::ReadMap("lattice.txt", "element.txt", parameter.map_filename_);
//  //}
//
//
//  auto supercell_config = Config::GenerateSupercell(parameter.supercell_size_,
//                                                    parameter.lattice_param_,
//                                                    "Al",
//                                                    parameter.structure_type_);
//
//
//  std::cout << "Finish config reading. Start CMC." << std::endl;
//
//  //auto cfg = Config::ReadPoscar("Al.poscar");
//  // auto supercell_cfg = Config::GenerateSupercell(4,4.046,"Al","FCC");
//  supercell_config.UpdateNeighborList({3.5, 4.8, 5.3});
//
//  auto atomVector = config.GetAtomVector();
//
//  std::set<Element> unique_set(atomVector.begin(), atomVector.end());
//  for (const auto& element : unique_set) {
//        std::cout << element.GetElementString() << " ";
//    }
//
//  PotentialEnergyEstimator
//  estimator("quartic_coefficients.json",config,supercell_config,
//            unique_set,3,3);
////
//  auto temp = estimator.GetEncodeVectorOfCluster(config,{1,4});
//  // std::cout << temp.size() << std::endl;
//
//  auto ele = atomVector[0];
//  std::cout <<  ele << std::endl;
//
//
//
//  std::cout << "element_set: ";
//  for (const auto& element : unique_set) {
//        std::cout << element.GetElementString() << " ";
//  }
//  std::cout << std::endl;
//
//}

// // #include "Config.h"
//
//
//
//
// int main(){
//
//   Config cfg = Config::ReadPoscar("Ti_Ni_BCC.poscar");
//   cfg.UpdateNeighborList({3.20,4.6,5.4}); // for BCC Ti_Ni
//
//   auto temp = cfg.GetNeighborAtomIdVectorOfAtom(4,1);
//
//   for (const auto& tmp : temp){
//     std::cout << tmp << std::endl;
//   }
//
//   size_t bond_order = 3;
//   cfg.WriteLattice("lattice_new.txt",bond_order);
//
//
//
// }
// // #include "Symmetry.h"
// int main() {

// Config cfg = Config::ReadPoscar("Al");
//  Config cfg = Config::ReadCfg("lowest_energy.cfg");
//  // FindPointGroupSymmetry(cfg);
//  // Config::WriteConfig("test.cfg", cfg);
//  // Config cfg = Config::ReadPoscar("Ti_Ni_BCC.poscar");
//
// cfg.UpdateNeighborList({3.5, 4.8, 5.3, 5.9, 6.5, 7.1, 7.6, 8.2});
// std::cout << "Neighbors updated" << std::endl;

// cfg.UpdateNeighborList({3.5, 4.8, 5.3, 5.9, 6.5, 7.1, 7.6, 8.2});
//  cfg.UpdateNeighborList({3.16, 4.47, 5.2, 5.26, 5.48, 5.92, 6.87, 7.87}); // for BCC tungsten

// Config cfg = Config::ReadPoscar("Ti_Ni_BCC.poscar");
// cfg.UpdateNeighborList({3.20,4.6,5.4}); // for BCC Ti_Ni
// PotentialEnergyEstimator
//     estimator("quartic_coefficients.json", cfg,
//               std::set<Element>{Element("Ti"),Element("Ni")},
//               3, 3);
// vector<size_t> lattice_id_pair = { 1, 2 };
// cout << "Configuration before swap :" << endl;
// auto start = estimator.GetEncodeVectorOfCluster(cfg, lattice_id_pair);
// std::pair<size_t,size_t> lattice_id_jump_pair = { 1, 2 };
// cfg.AtomJump(lattice_id_jump_pair);
// cout << "Configuration after swap :" << endl;
// auto end = estimator.GetEncodeVectorOfCluster(cfg, lattice_id_pair);

// PrintNeighborLists(cfg);

// PotentialEnergyEstimator
//     estimator("quartic_coefficients_rephrase.json", cfg,
//               std::set<Element>{Element("Mg"), Element("Al"), Element("Zn"), Element("Sn"), Element("X")},
//               3, 3);
// auto encode = estimator.GetEncodeVector(cfg);
// for (const auto &count : encode) {
//   std::cout << std::fixed << std::setprecision(0) << count << '\t';
// }
// std::cout << std::endl;
// std::cout << "The size is " << encode.size() << std::endl;

// cout << cfg.GetBasis() << endl;

// auto distance_order = cfg.GetDistanceOrder(1,3);

// PrintAtomVector(cfg);

// cout << "distance_order: " << distance_order << endl;

// auto tmp1 = InitializeLatticeClusterTypeSet(cfg, 3, 3);
// for (const auto &type : tmp1) {
//   std::cout << type << std::endl;
// }
//
// auto tmp3 = InitializeAtomClusterTypeSet(std::set<Element>{Element("Mg"), Element("Al"), Element("Zn")}, 3);
// for (const auto &type : tmp3) {
//   std::cout << type << std::endl;
// }

// auto tmp = InitializeClusterTypeSet(cfg,std::set<Element>{Element("W")},3,3);

// auto tmp2 = FindAllLatticeClusters(cfg ,3 ,3 ,{});
//
// for(const auto& cluster: tmp2){
//   cout << cluster.GetClusterType() << IdentifyAtomClusterType(cfg,cluster.GetLatticeIdVector()) << endl;
// }
//
// auto temp = estimator.GetEncodeVector(cfg);

//

//
// std:: cout  << "The Size of Cluster Type Set is " << tmp2.size() << endl;
//
// std::cout << std::endl;
// std::cout << "The Size of Encode Vector is " << encode.size() << std::endl;
//

// auto site1 = FindAllLatticeClusters2(cfg, 3, 3, {14191});

//

// cfg.GetNeighborLists();

// cout << endl;

// const auto &neighbor_lists = cfg.GetNeighborLists();
// std::vector<std::vector<size_t>> new_clusters;
// std::vector<std::vector<size_t>> old_clusters{{14191}};
//
//
//
//
// for (const auto &old_cluster : old_clusters) {
//  std::set<size_t> neighbors{};
//  std::cout << "Processing cluster: ";
//  for (auto id : old_cluster) {
//      std::cout << id << " ";
//  }
//  std::cout << std::endl;
//
//  // Retrieve neighbors up to max_bond_order
//  for (size_t m = 0; m < 3; m++) {
//      std::cout << "Bond order " << m + 1 << " neighbors: ";
//      for (auto lattice_id : old_cluster) {
//          const auto &current_neighbors = neighbor_lists[m][lattice_id];
//          for (auto neighbor_id : current_neighbors) {
//              neighbors.insert(neighbor_id);
//              std::cout << neighbor_id << cfg.GetElementOfLattice(neighbor_id) << " ";
//          }
//      }
//      std::cout << std::endl;
//  }
//
//  // Remove sites already in the cluster from neighbors
//  for (auto lattice_id : old_cluster) {
//      neighbors.erase(lattice_id);
//  }
//
//  std::cout << "Unique neighbors (after removal of current cluster sites): ";
//  for (auto neighbor : neighbors) {
//      std::cout << neighbor << " ";
//  }
//  std::cout << std::endl;
//
//  // Form new clusters by adding one neighbor to the current cluster
//  for (auto new_lattice_id : neighbors) {
//      std::vector<size_t> new_cluster{old_cluster};
//
//      // Check conditions for forming a valid cluster
//      if (std::all_of(old_cluster.begin(), old_cluster.end(), [&](size_t old_lattice_id) {
//          auto distance_order = cfg.GetDistanceOrder(old_lattice_id, new_lattice_id);
//          std::cout << "Checking distance order between " << old_lattice_id << " and "
//                    << new_lattice_id << ": " << distance_order << std::endl;
//          return old_lattice_id < new_lattice_id && distance_order <= 3;
//      })) {
//          new_cluster.push_back(new_lattice_id);
//          new_clusters.push_back(new_cluster);
//          std::cout << "New cluster formed: ";
//          for (auto id : new_cluster) {
//              std::cout << id << " ";
//          }
//          std::cout << std::endl;
//      }
//  }
//}
//
//

// cout << "Starting Configuration before swap for site2 :" << endl;
//
// auto start_site2 = estimator.GetEncodeVectorCluster(cfg, {13935});
//
// cout << endl;
//
//
//
//
//
// cout << endl;
//
//
// cout << "Starting Configuration after swap for site2 :" << endl;
//
// auto end_site2 = estimator.GetEncodeVectorCluster(cfg, {13935});
//
// cout << endl;

//

// auto cfg_supercell = Config::GenerateSupercell(4,4.046,"Al","FCC");
// cfg_supercell.UpdateNeighborList({3.5, 4.8, 5.3, 5.9, 6.5, 7.1, 7.6, 8.2});
//
// PotentialEnergyEstimator
//      estimator("quartic_coefficients.json", cfg,cfg_supercell,
//                std::set<Element>{Element("Al"),Element("Mg"),Element("Zn"),Element("X")},
//                3, 3);
//
// auto start_supercell = estimator.GetEncodeVector(cfg);
//
// cout << start_supercell << endl;

//
// auto temp = InitializeLatticeClusterTypeSet(cfg_supercell,3,3);
//
// auto temp2 = ConvertLatticeSetToHashMap(temp);
//
// for(const auto& type : temp2){
//   cout << type.first << " : " << type.second << endl;
// }
//
// auto temp3 = CountLatticeClusterTypes(cfg_supercell,3,3);
//
// for(const auto& type : temp3){
//   cout << type.first << ":" << type.second << endl;
// }

// std::cout << "Size of temp" << temp.size() << std::endl;

// return 0;

// }

// std::cout << "Distance Order: " << cfg.GetDistanceOrder(0, 1) << std::endl;
// Config cfg = Config::ReadPoscar("POSCAR30.gz");
// // cfg.SetPeriodicBoundaryCondition({false, false, false});

// std::vector<size_t> a(cfg.GetNumAtoms(), 0);
// std::vector<size_t> b(cfg.GetNumAtoms(), 0);
// std::vector<size_t> c(cfg.GetNumAtoms(), 0);
// std::vector<size_t> d(cfg.GetNumAtoms(), 0);
// std::vector<size_t> e(cfg.GetNumAtoms(), 0);
// std::vector<size_t> f(cfg.GetNumAtoms(), 0);
// std::vector<size_t> g(cfg.GetNumAtoms(), 0);
//
// for (size_t i = 0; i < cfg.GetNumAtoms(); ++i) {
//   a[i] = cfg.GetNeighborLists()[0][i].size();
//   b[i] = cfg.GetNeighborLists()[1][i].size();
//   c[i] = cfg.GetNeighborLists()[2][i].size();
//   d[i] = cfg.GetNeighborLists()[3][i].size();
//   e[i] = cfg.GetNeighborLists()[4][i].size();
//   f[i] = cfg.GetNeighborLists()[5][i].size();
//   g[i] = cfg.GetNeighborLists()[6][i].size();
// }
// std::map<std::string, Config::VectorVariant>
//     auxiliary_lists{
//     {"num_nn1", a},
//     {"num_nn2", b},
//     {"num_nn3", c},
//     {"num_nn4", d},
//     {"num_nn5", e},
//     {"num_nn6", f},
//     {"num_nn7", g}
// };
// cfg.WriteXyzExtended("output.xyz.gz", auxiliary_lists, {});
// cfg.WriteXyzExtended("output.xyz", auxiliary_lists, {});
// cfg.WriteConfig("output.cfg.gz");
// cfg.WriteConfig("output.cfg");

// // #include "Home.h"
// int main(int argc, char *argv[]) {
//   if (argc == 1) {
//     std::cout << "No input parameter filename." << std::endl;
//     return 1;
//   }
//   Parameter parameter(argc, argv);
//   Print(parameter);
//   Run(parameter);
//   return 0;
// }
