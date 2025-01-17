/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/18/23 4:41 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 9/27/23 11:29 AM                                                          *
 **************************************************************************************************/

#include "Symmetry.h"
#include "spglib.h"
#include <vector>
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <cmath>
#include <utility>
#include <spglib.h>
#include <iostream>
#include <vector>

void FindPointGroupSymmetry(const Config &config, std::vector<size_t> cluster) {

  SpglibDataset *dataset;

  double lattice_c[3][3];
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      lattice_c[i][j] = config.GetBasis()(static_cast<int>(i), static_cast<int>(j));



  size_t num_atom = cluster.size();
  double (*position)[3] = new double[num_atom][3];
  // Fill position array with atom coordinates
    for (size_t i = 0; i < num_atom; ++i) {
      std::cout << cluster[i] << std::endl;
        auto atom_position = config.GetRelativePositionOfLattice(cluster[i]);
        for (size_t j = 0; j < 3; ++j) {
            position[i][j] = atom_position[j];
        }
    }
  
  
  
  std::unordered_map<std::string, int> element_map = {{"Al", 1},{"Mg", 2},{"Zn", 3}};

  

  int *types = new int[num_atom];
  //for (size_t i = 0; i < num_atom; i++) { types[i] = 1; }
  for (size_t i = 0; i < num_atom; i++){ 
    auto element = config.GetElementOfLattice(cluster[i]);

    types[i] = element_map[element.GetElementString()];
  }
   
  
  // SpglibDataset has to be freed after use.
  dataset = spg_get_dataset(lattice_c, position, types, static_cast<int>(num_atom), 1e-5);
  // delete[] position;
  // delete[] types;

  // printf("International symbol: %s (%d)\n", dataset->international_symbol,
  //        dataset->spacegroup_number);
  // printf("Hall symbol:   %s\n", dataset->hall_symbol);
  // printf("Wyckoff letters:\n");
// 
  // printf("\n");
  printf("Equivalent atoms:\n");
  // for (size_t i = 0; i < dataset->n_atoms; i++) {
  //   printf("%zu -> %d\n", i, dataset->equivalent_atoms[i]);
  // }


  std::vector<int> eq_vector;

  for(size_t i=0; i<dataset->n_atoms; i++){
    std::cout << i << "; " << cluster[i] << "; " << dataset->equivalent_atoms[i] << std::endl;
    eq_vector.push_back(dataset->equivalent_atoms[i]);
  }

  std::cout << "Equivalent Atom Vector: {";
  for(auto id : eq_vector){
    std::cout << id << ", ";
  }
  std::cout << "}"  << std::endl;


  std::set<int> eq_set(eq_vector.begin(), eq_vector.end());

  std::vector<std::vector<int>> temp_vector;

  for (auto idx : eq_set) {
    std::vector<int> temp_v;
    for (int i = idx; i<eq_vector.size()+1; ++i) {
      if (idx != eq_vector[i]) {
        temp_vector.push_back(temp_v);
        break;
      }
      else{
        temp_v.push_back(i);
        std::cout << i << std::endl;
      }
    }

  }


  std::cout << "MMM vector: {";
  for (auto vec : temp_vector) {
    std::cout << "{ ";
    for (auto i : vec) {
      std::cout << i << ", ";
    }
    std::cout << "}, ";
  }
  std::cout << "}" << std::endl;

  


  printf("Point group symbol: %s\n", dataset->pointgroup_symbol);

  std::cout << "Lattice ID: {";

  for(const auto& id : cluster){
    std::cout << id << config.GetElementOfLattice(id) << ", ";
  }
  std::cout << "} ; ";

  if (cluster.size() == 2){
    std::cout << "Distance Order: " << config.GetDistanceOrder(cluster[0],cluster[1]) << "; ";
  }
  std::cout << "Point Group Symbol: " << dataset->pointgroup_symbol << std::endl;

  // Deallocate SpglibDataset, otherwise induce memory leak.
  spg_free_dataset(dataset);
}



void FindPointGroupSymmetryMM2(const Config &config, std::vector<size_t> cluster) {

  SpglibDataset *dataset;

  double lattice_c[3][3];
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      lattice_c[i][j] = config.GetBasis()(static_cast<int>(i), static_cast<int>(j));



  size_t num_atom = cluster.size();
  double (*position)[3] = new double[num_atom][3];
  // Fill position array with atom coordinates
    for (size_t i = 0; i < num_atom; ++i) {
      std::cout << cluster[i] << std::endl;
        Eigen::Vector3d atom_position = config.GetRelativePositionOfLattice(cluster[i]);
        for (size_t j = 0; j < 3; ++j) {
            position[i][j] = atom_position[j];
        }
    }



  position[16][1] += 0.01;
  
  // position[0][2] += 0.01;
  std::unordered_map<std::string, int> element_map = {{"Al", 1},{"Mg", 2},{"Zn", 3}};

  

  int *types = new int[num_atom];
  //for (size_t i = 0; i < num_atom; i++) { types[i] = 1; }
  for (size_t i = 0; i < num_atom; i++){ 
    auto element = config.GetElementOfLattice(cluster[i]);

    types[i] = element_map[element.GetElementString()];
  }
   
  
  // SpglibDataset has to be freed after use.
  dataset = spg_get_dataset(lattice_c, position, types, static_cast<int>(num_atom), 1e-5);
  // delete[] position;
  // delete[] types;

  // printf("International symbol: %s (%d)\n", dataset->international_symbol,
  //        dataset->spacegroup_number);
  // printf("Hall symbol:   %s\n", dataset->hall_symbol);
  // printf("Wyckoff letters:\n");
// 
  // printf("\n");
  printf("Equivalent atoms:\n");
  // for (size_t i = 0; i < dataset->n_atoms; i++) {
  //   printf("%zu -> %d\n", i, dataset->equivalent_atoms[i]);
  // }


  std::vector<int> eq_vector;

  for(size_t i=0; i<dataset->n_atoms; i++){
    std::cout << i << "; " << cluster[i] << "; " << dataset->equivalent_atoms[i] << std::endl;
    eq_vector.push_back(dataset->equivalent_atoms[i]);
  }

  std::cout << "Equivalent Atom Vector: {";
  for(auto id : eq_vector){
    std::cout << id << ", ";
  }
  std::cout << "}"  << std::endl;


  std::set<int> eq_set(eq_vector.begin(), eq_vector.end());

  std::vector<std::vector<int>> temp_vector;

  for (auto idx : eq_set) {
    std::vector<int> temp_v;
    for (int i = idx; i<eq_vector.size()+1; ++i) {
      if (idx != eq_vector[i]) {
        temp_vector.push_back(temp_v);
        break;
      }
      else{
        temp_v.push_back(i);
        std::cout << i << std::endl;
      }
    }

  }


  std::cout << "MM2 vector: {";
  for (auto vec : temp_vector) {
    std::cout << "{ ";
    for (auto i : vec) {
      std::cout << i << ", ";
    }
    std::cout << "}, ";
  }
  std::cout << "}" << std::endl;

  


  printf("Point group symbol: %s\n", dataset->pointgroup_symbol);

  std::cout << "Lattice ID: {";

  for(const auto& id : cluster){
    std::cout << id << config.GetElementOfLattice(id) << ", ";
  }
  std::cout << "} ; ";

  if (cluster.size() == 2){
    std::cout << "Distance Order: " << config.GetDistanceOrder(cluster[0],cluster[1]) << "; ";
  }
  std::cout << "Point Group Symbol: " << dataset->pointgroup_symbol << std::endl;

  // Deallocate SpglibDataset, otherwise induce memory leak.
  spg_free_dataset(dataset);
}



