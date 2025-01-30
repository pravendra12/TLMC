#include "VacancyMigrationPredictor.h"
#include <fstream>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

std::unordered_map<std::string, double>
ReadParamsFromJson(const std::string &json_filename) {
  std::ifstream ifs(json_filename, std::ifstream::in);
  if (!ifs.is_open()) {
        throw std::runtime_error("Failed to open JSON file.");
    }
  json all_parameters;
  ifs >> all_parameters;

  std::unordered_map<std::string, double> param_map;

  if (all_parameters.contains("parameters")) {
    auto parameters = all_parameters["parameters"];

    for (auto& [key, value] : parameters.items()) {
        param_map[key] = value;
    }
  } 
  else {
    throw  std::runtime_error("No 'parameters' are present in the JSON file");
  }

  return param_map;
}


VacancyMigrationPredictor::VacancyMigrationPredictor(
          const std::string &predictor_filename){
  
  auto parameter_map = ReadParamsFromJson(predictor_filename);
  // Reading parameters from predictor file.
  try {
    slope_E_diff_ = parameter_map.at("slope_E_diff");
    intercept_E_diff_ = parameter_map.at("intercept_E_diff");
    slope_DbE_ = parameter_map.at("slope_DbE");
    slope_DbA_ = parameter_map.at("slope_DbA");
    intercept_barrier_ = parameter_map.at("intercept_barrier");
  } catch (...){
    throw std::runtime_error("All parameters are not present for Vacancy Migration Predictor");
  }
}


// std::pair<double, double> 
// VacancyMigrationPredictor::GetBarrierAndDiffFromLatticeIdPair(
//                  Config &config,
//                  const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
// 
//   auto dE = GetDiff(config, lattice_id_jump_pair);
//   auto barrier = GetBarrier(config, lattice_id_jump_pair);
//   std::pair<double, double> barrier_de = {barrier,dE};
//   
//   return barrier_de;
// }

std::array<double, 3>
VacancyMigrationPredictor::GetBarrierAndDiffFromLatticeIdPair(
                 Config &config,
                 const std::pair<size_t, size_t> &lattice_id_jump_pair) const{
  
  std::array<double, 3> barrier_de;

  // Element at first lattice ID moves to second lattice ID

  auto forward_barrier = GetBarrierNew(config, lattice_id_jump_pair);
  auto forward_Ed = GetDiffNew(config, lattice_id_jump_pair);
  
  std::pair<size_t, size_t> backward_jump_pair = {lattice_id_jump_pair.second, 
                                                  lattice_id_jump_pair.first};

  auto backward_barrier = GetBarrierNew(config, backward_jump_pair);
  auto backward_Ed = GetDiffNew(config, backward_jump_pair);

  double new_forward_barrier;
  double new_forward_Ed;
  double new_backward_barrier;
  double new_backward_Ed;

  double average_barrier;
  double average_Ed = (abs(forward_Ed) + abs(backward_Ed))/2;
  
  if (forward_barrier < backward_barrier) {
    average_barrier = (abs(forward_barrier) + 
                      (abs(backward_barrier)-abs(backward_Ed)))/2;

    new_forward_barrier = average_barrier;
    new_forward_Ed = -average_Ed;
    new_backward_barrier = average_barrier + average_Ed;
    new_backward_Ed = average_Ed;
  }
  else {
    average_barrier = ((abs(forward_barrier) - abs(forward_Ed)) + 
                        abs(backward_barrier))/2;

    new_forward_barrier = average_barrier + average_Ed;
    new_forward_Ed = average_Ed;
    new_backward_barrier = average_barrier ;
    new_backward_Ed = -average_Ed;
  }

  barrier_de[0] = new_forward_barrier;
  barrier_de[1] = new_backward_barrier;
  barrier_de[2] = new_forward_Ed;

  if (new_backward_Ed != -new_forward_Ed) {
    std::cout << "Hey Man there is some error please look into it" << std::endl;
  }
  
  return barrier_de;
}


double VacancyMigrationPredictor::GetDiff(Config &config, 
                  const std::pair<size_t, size_t> lattice_id_jump_pair) const {
  
  // According to the paper they have emphasised about the vacancy migration
  // Whereas in CMC we swap two atoms randomly, not sure whether it will work 
  // in that case. Therefore, for this function is specifically for computing 
  // driving force for vacancy migration.
  
  size_t migrating_lattice_id; // final vacancy position
  size_t vacancy_id; // initial vacancy position 
  
  // Getting the lattice ID of migration atom and vacancy
  if (config.GetVacancyLatticeId() == lattice_id_jump_pair.first) {
    vacancy_id = lattice_id_jump_pair.first;
    migrating_lattice_id = lattice_id_jump_pair.second;
  }
  else {
    vacancy_id = lattice_id_jump_pair.second;
    migrating_lattice_id = lattice_id_jump_pair.first;
  }
 
  // Df for initial vacancy structure
  auto Df_initial = GetDf(config, vacancy_id);

  // Df for final vacancy structure
  
  auto migrating_lattice_id_neighbours = config.GetNeighborLatticeIdVectorOfLattice(migrating_lattice_id, 1);
  auto migrating_lattice_element = config.GetElementOfLattice(migrating_lattice_id);

  std::vector<double> xSv_vector;
  for (auto &id : migrating_lattice_id_neighbours) {
    auto element = config.GetElementOfLattice(id);
    if (element.GetElementString() == "X") {
      xSv_vector.push_back(GetxSv(migrating_lattice_element));
  //    std::cout << "Put the Element in place of the Vacancy" << std::endl;
    }
    else {
      xSv_vector.push_back(GetxSv(element));
    }
  }
//
  auto Df_final = GetGeometricMean(xSv_vector, 1);

//  config.LatticeJump({vacancy_id, migrating_lattice_id});

  // Df at final vacancy position, which is atom's initial position
  // auto Df_final = GetDf(config, migrating_lattice_id);

//  config.LatticeJump({vacancy_id, migrating_lattice_id});
  
  // std::cout << vacancy_id << config.GetElementOfLattice(vacancy_id) << " : " <<
  // atom_id << config.GetElementOfLattice(atom_id) << std::endl;

  return slope_E_diff_*(Df_final - Df_initial) + intercept_E_diff_;

}

double VacancyMigrationPredictor::GetBarrier(
                const Config &config,
                const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  
  // Assuming {k_lattice_id, i_lattice_id}
  // k_lattice_id -> Vacancy Lattice Id
  // i_lattice_id -> Migration Atoms Lattice Id
  // Below method takes care of that

  size_t atom_id;
  size_t vacancy_id;
  
  // Getting the lattice ID of migration atom and vacancy
  if (config.GetVacancyLatticeId() == lattice_id_jump_pair.first) {
    vacancy_id = lattice_id_jump_pair.first;
    atom_id = lattice_id_jump_pair.second;
  }
  else {
    vacancy_id = lattice_id_jump_pair.second;
    atom_id = lattice_id_jump_pair.first;
  }

  // Classification of 1st NN Atoms as I1, I2, F1, F2 based on the order

  size_t bond_order = 1; // Only considering 1st nearest neighbours
  
  auto neighboring_lattice_ids = config.GetNeighboringLatticeIdSetOfPair(
                                            lattice_id_jump_pair, bond_order);

  std::vector<double> i1_positions; // I1
  std::vector<double> i2_positions; // I2
  std::vector<double> f1_positions; // F1
  std::vector<double> f2_positions; // F2

  for (auto& id : neighboring_lattice_ids) {
    // bond order between id and migration atom
    auto initial_BO = config.GetDistanceOrder(atom_id, id);
    // bond order between id and vacancy, which will be final position of the atom
    auto final_BO = config.GetDistanceOrder(vacancy_id, id);

    auto atom_element = config.GetElementOfLattice(id);

    // condition for I1 atoms
    if (initial_BO == 1 && final_BO == 2) {
      i1_positions.push_back(GetxSv(atom_element));
    }    

    // condition for F1 atoms
    else if (initial_BO == 2 && final_BO == 1) {
      f1_positions.push_back(GetxSv(atom_element));
    }

    // condition for I2 atoms
    else if (initial_BO == 1 && (final_BO >= 3)) {
      i2_positions.push_back(GetxSv(atom_element));
    }

    // condition for F2 atoms
    else if ((initial_BO >= 3) && final_BO == 1 ) {
      f2_positions.push_back(GetxSv(atom_element));
    }
    
  }
  
  // Effect of Migrating Atom
  double DbA = GetxSv(config.GetElementOfLattice(atom_id)); 

  auto i1_xSv = GetGeometricMean(i1_positions, 1.0);
  auto i2_xSv = GetGeometricMean(i2_positions, 1.5);
  auto f1_xSv = GetGeometricMean(f1_positions, -1.5);
  
  // Effect of Migrating Environment
  double DbE =  i1_xSv*i2_xSv*f1_xSv; 
  
  return slope_DbE_*DbE + slope_DbA_*DbA + intercept_barrier_;
}

double GetGeometricMean(std::vector<double> &xSv_vector, double z = 1.0) {
  double product = 1;
  double n = static_cast<double>(xSv_vector.size());

  for (auto& xSv : xSv_vector) {
    product *= xSv;
  }

  return pow(product, z/n);
}

double GetxSv(const Element &element) {
  return element.GetElectronegativity() * element.GetSv();
}

double GetDf(const Config &config, const size_t lattice_id) {
  // Only first NN
  auto neighbouring_lattice_ids = config.GetNeighborLatticeIdVectorOfLattice(lattice_id, 1);

  std::vector<double> xSv_vector{};

  for (auto &id : neighbouring_lattice_ids) {
    auto element = config.GetElementOfLattice(id);
    xSv_vector.push_back(GetxSv(element));
  }

  return GetGeometricMean(xSv_vector);
}




double VacancyMigrationPredictor::GetBarrierNew(
                const Config &config,
                const std::pair<size_t, size_t> &lattice_id_jump_pair) const {
  
  // Assuming {k_lattice_id, i_lattice_id}
  // k_lattice_id -> Vacancy Lattice Id
  // i_lattice_id -> Migration Atoms Lattice Id
  // Below method takes care of that
  
  // position of atom
  size_t initial_position = lattice_id_jump_pair.first;
  size_t final_position = lattice_id_jump_pair.second;

  // Classification of 1st NN Atoms as I1, I2, F1, F2 based on the order

  size_t bond_order = 1; // Only considering 1st nearest neighbours
  
  auto neighboring_lattice_ids = config.GetNeighboringLatticeIdSetOfPair(
                                            lattice_id_jump_pair, bond_order);

  std::vector<double> i1_positions; // I1
  std::vector<double> i2_positions; // I2
  std::vector<double> f1_positions; // F1
  std::vector<double> f2_positions; // F2

  for (auto& id : neighboring_lattice_ids) {
    // bond order between id and migration atom
    auto initial_BO = config.GetDistanceOrder(initial_position, id);
    // bond order between id and vacancy, which will be final position of the atom
    auto final_BO = config.GetDistanceOrder(final_position, id);

    auto atom_element = config.GetElementOfLattice(id);

    // condition for I1 atoms
    if (initial_BO == 1 && final_BO == 2) {
      i1_positions.push_back(GetxSv(atom_element));
    }    

    // condition for F1 atoms
    else if (initial_BO == 2 && final_BO == 1) {
      f1_positions.push_back(GetxSv(atom_element));
    }

    // condition for I2 atoms
    else if (initial_BO == 1 && (final_BO >= 3)) {
      i2_positions.push_back(GetxSv(atom_element));
    }

    // condition for F2 atoms
    else if ((initial_BO >= 3) && final_BO == 1 ) {
      f2_positions.push_back(GetxSv(atom_element));
    }
    
  }

  double DbA;
  
  // Effect of Migrating Atom
  auto element_at_initial_position = config.GetElementOfLattice(initial_position);
  auto element_at_final_position = config.GetElementOfLattice(final_position);
  
  if (element_at_initial_position.GetElementString() == "X") {
    DbA = GetxSv(element_at_final_position); 
  }
  else {
    DbA = GetxSv(element_at_initial_position);
  }


  auto i1_xSv = GetGeometricMean(i1_positions, 1.0);
  auto i2_xSv = GetGeometricMean(i2_positions, 1.5);
  auto f1_xSv = GetGeometricMean(f1_positions, -1.5);
  
  // Effect of Migrating Environment
  double DbE =  i1_xSv*i2_xSv*f1_xSv; 
  
  return slope_DbE_*DbE + slope_DbA_*DbA + intercept_barrier_;
}



double VacancyMigrationPredictor::GetDiffNew(Config &config, 
                        const std::pair<size_t, size_t> lattice_id_jump_pair) const {
  
  // According to the paper they have emphasised about the vacancy migration
  // Whereas in CMC we swap two atoms randomly, not sure whether it will work 
  // in that case. Therefore, for this function is specifically for computing 
  // driving force for vacancy migration.

  // Assumption is that atom is at first lattice Id and 
  // moves to second lattice Id where vacancy is there.

  double Df_initial;
  double Df_final;

  auto first_element = config.GetElementOfLattice(lattice_id_jump_pair.first);
  auto second_element = config.GetElementOfLattice(lattice_id_jump_pair.second);
  
  // When vacancy is at second lattice Id

  if (second_element.GetElementString() == "X") {
    
    Df_initial = GetDf(config, lattice_id_jump_pair.second);
    
    auto migrating_lattice_id_neighbours = 
      config.GetNeighborLatticeIdVectorOfLattice(lattice_id_jump_pair.first, 1);
    auto migrating_lattice_element = 
      config.GetElementOfLattice(lattice_id_jump_pair.first);
  
    std::vector<double> xSv_vector;
    for (auto &id : migrating_lattice_id_neighbours) {
      auto element = config.GetElementOfLattice(id);
      if (element.GetElementString() == "X") {
        xSv_vector.push_back(GetxSv(migrating_lattice_element));
      }
      else {
        xSv_vector.push_back(GetxSv(element));
      }
    }
    Df_final = GetGeometricMean(xSv_vector, 1);
  
  }
  
  // Now assuming the swap have happened but in reality 
  // vacancy is still at second position and atom at first position
  // To compute the driving force such that when vacancy is at first position (assumption)
  // And atom is at second.

  else {

    auto migrating_lattice_id_neighbours = 
    config.GetNeighborLatticeIdVectorOfLattice(lattice_id_jump_pair.second, 1);
    auto migrating_lattice_element = 
    config.GetElementOfLattice(lattice_id_jump_pair.second);

    std::vector<double> xSv_vector;
    for (auto &id : migrating_lattice_id_neighbours) {
      auto element = config.GetElementOfLattice(id);
      if (element.GetElementString() == "X") {
        xSv_vector.push_back(GetxSv(migrating_lattice_element)); 
      }
      else {
        xSv_vector.push_back(GetxSv(element));
      }
    }

    Df_initial = GetGeometricMean(xSv_vector, 1);
    
    Df_final = GetDf(config, lattice_id_jump_pair.first);
  }
  
  return slope_E_diff_*(Df_final - Df_initial) + intercept_E_diff_;

}
