/**************************************************************************************************
 * Copyright (c) 2020-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 1/16/20 3:55 AM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 10/25/23 10:35 PM                                                         *
 **************************************************************************************************/

/*! \file  main.cpp
 *  \brief File for the main function.
 */


// #include "ClusterExpansion.h"
// #include "PotentialEnergyEstimator.h"
// using namespace std;
// 
// // #include <chrono>
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

#include "Home.h"
// #include "ThermodynamicAveraging.h"
// #include "CanonicalMcSerial.h"
// #include "CanonicalMcAbstract.h"
// #include "PotentialEnergyEstimator.h"
// #include "Config.h"
// #include "ClusterExpansion.h"
// #include <unordered_set>
// #include <unordered_map>
// #include <vector>
// #include <string>
// #include <set>
// #include <algorithm>
// #include <boost/functional/hash.hpp>
// #include "CanonicalMcSerial.h"
// #include "Symmetry.h"
// #include <omp.h>
// #include <chrono>
// #include <Eigen/Dense>
// #include <vector>
// #include <algorithm>
#include "VacancyMigrationPredictor.h"
#include "JumpEvent.h"
// // #include "KineticMcChainOmpi.h"
// // #include "KineticMcAbstract.h"
// #include "JumpEvent.h"
// #include "Constants.hpp"
// #include "ShortRangeOrder.h"
// #include "Traverse.h"
// #include "KineticMcFirstMpi.h"
// #include "SymmetryCustom.h"
// #include "LatticeClusterMMM.hpp"
// #include "Element.hpp"
// #include <cmath>




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

int main(int argc, char *argv[]) {
  if (argc == 1) {
    std::cout << "No input parameter filename." << std::endl;
    return 1;
  }
  api::Parameter parameter(argc, argv);
  api::Print(parameter);
  api::Run(parameter);
}


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
//     // so this GetBarrierNew function takes the lattice jump pair and will return 
//     // barrier for the event where the first lattice ID wants to move to second 
//     // lattice ID ; one of them need to be vacancy
//     
//     // forward barrier will be the case where the element will move 
//     // from migrating atom id to the vacancy id
//     auto forward_barrier = migrationPredictor.GetBarrierNew(cfg, {migratingAtomId, vacancyId});
// 
//     auto forward_Ed = migrationPredictor.GetDiffNew(cfg, {migratingAtomId, vacancyId});
// 
//     std::cout << migratingAtomId << cfg.GetElementOfLattice(migratingAtomId) << " ";
//     std::cout << forward_barrier << " " << forward_Ed << std::endl;
// 
//     
//     std::cout << "Comparsion between the previous barrier and current " << 
//     "barrier func for atom moving to vacany position" << std::endl;
//     std::cout << "Old method barrier : " << migrationPredictor.GetBarrier(cfg, {migratingAtomId, vacancyId}) << std::endl;
//     std::cout << "New method barrier : " << forward_barrier  << std::endl;
// 
//     std::cout << "Old method Ed : " << migrationPredictor.GetDiff(cfg, {migratingAtomId, vacancyId}) << std::endl;
//     std::cout << "New method Ed : " << forward_Ed  << std::endl;
// 
// 
//     
// 
//     // how to get the barrier and Ed without swapping the positions
//       
//     // backward barrier will be case where your atom will be at the vacancy id 
//     // and will want to move back to the migrating atom id
//     auto backward_barrier = migrationPredictor.GetBarrierNew(cfg, {vacancyId, migratingAtomId});
//     auto backward_Ed = migrationPredictor.GetDiffNew(cfg, {vacancyId, migratingAtomId});
//     
//     cfg.LatticeJump({vacancyId, migratingAtomId});
//     
//     std::cout << "Comparsion between the backward previous barrier and current " << 
//     "barrier func for atom moving to vacany position" << std::endl;
// 
//     std::cout << "Old method barrier : " << migrationPredictor.GetBarrier(cfg, {migratingAtomId, vacancyId}) << std::endl;
//     std::cout << "New method barrier : " << backward_barrier  << std::endl;
// 
//     std::cout << "Old method Ed : " << migrationPredictor.GetDiff(cfg, {migratingAtomId, vacancyId}) << std::endl;
//     std::cout << "New method Ed : " << backward_Ed  << std::endl;
// 
// 
//     cfg.LatticeJump({vacancyId, migratingAtomId});
//     
//     std::cout << migratingAtomId << cfg.GetElementOfLattice(migratingAtomId) << " ";
//     std::cout << backward_barrier << " " << backward_Ed << std::endl;
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
//     std::cout << "Forward Barrier #: " << forward_backward_info.first.first << std::endl;
//     std::cout << "Forward Ed #: " << forward_backward_info.first.second << std::endl;
//     std::cout << "Backward Barrier #: " << forward_backward_info.second.first << std::endl;
//     std::cout << "Backward Ed #: " << forward_backward_info.second.second << std::endl;
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
//     std::cout << "Forward Barrier #: " << event.GetForwardBarrier() << std::endl;
//     std::cout << "Forward Ed #: " << event.GetEnergyChange() << std::endl;
//     std::cout << "Backward Barrier: " << event.GetBackwardBarrier() << std::endl;
// 
//     auto backward_event = event.GetReverseJumpEvent();
//     
//     std::cout << "forward Barrier  for b#: " << backward_event.GetForwardBarrier() << std::endl;
//     std::cout << " Ed for backward event #: " << backward_event.GetEnergyChange() << std::endl;
//     std::cout << "backward Barrier  for b#: " << backward_event.GetBackwardBarrier() << std::endl;
// 
//     break;
// 
//   }
//   
// }

//int main(int argc, char *argv[]) {
//  // if (argc == 1) {
//  //   std::cout << "No input parameter filename." << std::endl;
//  //   return 1;
//  // }
//  // api::Parameter parameter(argc, argv);
//  // api::Print(parameter);
//  // api::Run(parameter);
//
//  auto cfg = Config::ReadConfig("TiTaMoNb_Vac.POSCAR");
//  cfg.UpdateNeighborList({3.20, 4.6, 5.4});
//
//
//  auto old_vacancyId = cfg.GetVacancyLatticeId();
//
//  cfg.LatticeJump({old_vacancyId, cfg.GetNeighborLatticeIdVectorOfLattice(old_vacancyId, 1)[0]});
//  
//  auto vacancyId = cfg.GetVacancyLatticeId();
//
//
//  for (int i=0; i<8; i++){
//
//    auto migratingAtomId = cfg.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1)[i];
//
//    std::cout << "Lattice ID Pair : " << vacancyId << " " << migratingAtomId << std::endl;
//
//    VacancyMigrationPredictor migrationPredictor("predictor_file.json");
//
//    auto forward_barrier = migrationPredictor.GetBarrier(cfg, {vacancyId, migratingAtomId});
//    auto forward_Ed = migrationPredictor.GetDiff(cfg, {vacancyId, migratingAtomId});
//
//    std::cout << migratingAtomId << cfg.GetElementOfLattice(migratingAtomId) << " ";
//    std::cout << forward_barrier << " " << forward_Ed << std::endl;
//
//    cfg.LatticeJump({vacancyId, migratingAtomId});
//
//    // how to get the barrier and Ed without swapping the positions
//      
//
//    auto backward_barrier = migrationPredictor.GetBarrier(cfg, {vacancyId, migratingAtomId});
//    auto backward_Ed = migrationPredictor.GetDiff(cfg, {vacancyId, migratingAtomId});
//    
//    std::cout << migratingAtomId << cfg.GetElementOfLattice(migratingAtomId) << " ";
//    std::cout << backward_barrier << " " << backward_Ed << std::endl;
//
//    cfg.LatticeJump({vacancyId, migratingAtomId});
//
//
//
//
//    double new_forward_barrier;
//    double new_forward_Ed;
//    double new_backward_barrier;
//    double new_backward_Ed;
//
//    double average_barrier;
//    double average_Ed;
//
//    if (forward_barrier < backward_barrier) {
//
//      average_Ed = (abs(forward_Ed) + abs(backward_Ed))/2;
//      average_barrier = (abs(forward_barrier) + (abs(backward_barrier)-abs(backward_Ed)))/2;
//
//      std::cout << "Average Ed : " << average_Ed << std::endl;
//      std::cout << "Average barrier : " << average_barrier << std::endl; 
//      new_forward_barrier = average_barrier;
//      new_forward_Ed = -average_Ed;
//      new_backward_barrier = average_barrier + average_Ed;
//      new_backward_Ed = average_Ed;
//    }
//    else {
//      average_Ed = (abs(forward_Ed) + abs(backward_Ed))/2;
//      average_barrier = ((abs(forward_barrier) - abs(forward_Ed)) + abs(backward_barrier))/2;
//
//      std::cout << "Average Ed : " << average_Ed << std::endl;
//      std::cout << "Average barrier : " << average_barrier << std::endl; 
//      new_forward_barrier = average_barrier + average_Ed;
//      new_forward_Ed = average_Ed;
//      new_backward_barrier = average_barrier ;
//      new_backward_Ed = -average_Ed;
//
//    }
//
//    std::cout << "New barrier and Ed from the proposed solution" << std::endl;
//
//    std::cout << "Forward Barrier : " << new_forward_barrier << std::endl;
//    std::cout << "Forward Ed : " << new_forward_Ed << std::endl;
//    std::cout << "Backward Barrier : " << new_backward_barrier << std::endl;
//    std::cout << "Backward Ed : " << new_backward_Ed << std::endl;
//
//    std::cout << "-----------------------" << std::endl;
//
//    std::pair<std::pair<double, double>, std::pair<double, double>> forward_backward_info =
//    {{new_forward_barrier, new_forward_Ed}, {new_backward_barrier, new_backward_Ed}};
//
//    
//
//    std::cout << "Forward Barrier #: " << forward_backward_info.first.first << std::endl;
//    std::cout << "Forward Ed #: " << forward_backward_info.first.second << std::endl;
//    std::cout << "Backward Barrier #: " << forward_backward_info.second.first << std::endl;
//    std::cout << "Backward Ed #: " << forward_backward_info.second.second << std::endl;
//
//    std::array<double, 3> barrier_de;
//    barrier_de[0] = new_forward_barrier;
//    barrier_de[1] = new_backward_barrier;
//    barrier_de[2] = new_forward_Ed;
//
//    double beta = 1 / 400 / constants::kBoltzmann;
//
//    mc::JumpEvent event({vacancyId, migratingAtomId},
//                        barrier_de,
//                        beta);
//    
//    std::cout << "Defining the event " << std::endl;
//    std::cout << std::endl;
//    
//    std::cout << "For Forward event " << std::endl;
//
//    std::cout << "Forward Barrier #: " << event.GetForwardBarrier() << std::endl;
//    std::cout << "Forward Ed #: " << event.GetEnergyChange() << std::endl;
//    std::cout << "Backward Barrier: " << event.GetBackwardBarrier() << std::endl;
//
//    auto backward_event = event.GetReverseJumpEvent();
//    
//    std::cout << "forward Barrier  for b#: " << backward_event.GetForwardBarrier() << std::endl;
//    std::cout << " Ed for backward event #: " << backward_event.GetEnergyChange() << std::endl;
//    std::cout << "backward Barrier  for b#: " << backward_event.GetBackwardBarrier() << std::endl;
//
//  }
//  
//}
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



    //auto sorted_ids = SymmetricallySortLatticeIDs(cfg,temp_vector);





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
    

    //FindEquivalentClusters(cfg,pair_cluster_vector);
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
  
  
  
  //auto rot_mat = GetLatticePairRotationMatrix(cfg,lattice_id_pair);
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
  
    //auto end1 = std::chrono::high_resolution_clock::now();
//
    //auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start).count();
//
    //std::cout << "Time taken for mc_simulation to initialize: " << duration1 << " milliseconds" << std::endl;

    // mc_simulation.Simulate();
    

    
    // auto temp = IndentifyLatticeClusterType(cfg,{1379,1650,1072});
    // std::cout << temp << std::endl;
    // cfg.GetDistanceOrder(1379,1650);
    // cfg.GetDistanceOrder(1072,1650);
    // cfg.GetDistanceOrder(1379,1072);



    //printNestedList(n_list_lattice);

    //std::cout << "Distance Order: " << cfg.GetDistanceOrder(1090,1667) << std::endl;
    
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
//auto nn_711 = cfg.GetNeighborLatticeIdVectorOfLattice(711, 1);
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

// #include "Config.h"
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
// #include "Symmetry.h"
// int main() {

  //Config cfg = Config::ReadPoscar("Al");
  // Config cfg = Config::ReadCfg("lowest_energy.cfg");
  // // FindPointGroupSymmetry(cfg);
  // // Config::WriteConfig("test.cfg", cfg);
  // // Config cfg = Config::ReadPoscar("Ti_Ni_BCC.poscar");
// 
  // cfg.UpdateNeighborList({3.5, 4.8, 5.3, 5.9, 6.5, 7.1, 7.6, 8.2});
  // std::cout << "Neighbors updated" << std::endl;

  //cfg.UpdateNeighborList({3.5, 4.8, 5.3, 5.9, 6.5, 7.1, 7.6, 8.2});
  // cfg.UpdateNeighborList({3.16, 4.47, 5.2, 5.26, 5.48, 5.92, 6.87, 7.87}); // for BCC tungsten
  
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

  //auto tmp = InitializeClusterTypeSet(cfg,std::set<Element>{Element("W")},3,3);

 

  

  

  // auto tmp2 = FindAllLatticeClusters(cfg ,3 ,3 ,{});
// 
  // for(const auto& cluster: tmp2){
  //   cout << cluster.GetClusterType() << IndentifyAtomClusterType(cfg,cluster.GetLatticeIdVector()) << endl;
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
  


  //auto site1 = FindAllLatticeClusters2(cfg, 3, 3, {14191});
  
  

  //  

  //cfg.GetNeighborLists();

  // cout << endl;
  
  //const auto &neighbor_lists = cfg.GetNeighborLists();
  //std::vector<std::vector<size_t>> new_clusters;
  //std::vector<std::vector<size_t>> old_clusters{{14191}};
//
//
//
//
  //for (const auto &old_cluster : old_clusters) {
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

// #include "Home.h"
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
