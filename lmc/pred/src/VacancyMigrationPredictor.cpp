#include "VacancyMigrationPredictor.h"



VacancyMigrationPredictor::VacancyMigrationPredictor(const std::string &predictor_filename) {
  
  // To be done in future, slope and intercepts to predict the barrier and
  // driving force will be read from the predictor file, for now using some,
  // random values.

}

std::pair<double, double> 
VacancyMigrationPredictor::GetBarrierAndDiffFromLatticeIdPair(
                 const Config &config,
                 const std::pair<size_t, size_t> &lattice_id_jump_pair) const {

  auto dE = GetDiff(config, lattice_id_jump_pair);
  auto barrier = GetBarrier(config, lattice_id_jump_pair);
  std::pair<double, double> barrier_de = {barrier,dE};
  
  return barrier_de;
}

double VacancyMigrationPredictor::GetDiff(const Config &config, 
                        const std::pair<size_t, size_t> lattice_id_jump_pair) const {
  
  // Assuming atom at first lattice Id wants to move to second lattice Id
  
  size_t atom_id; // final vacancy position
  size_t vacancy_id; // initial vacancy position 
  
  // Getting the lattice ID of migration atom and vacancy
  if (config.GetVacancyLatticeId() == lattice_id_jump_pair.first) {
    vacancy_id = lattice_id_jump_pair.first;
    atom_id = lattice_id_jump_pair.second;
  }
  else {
    vacancy_id = lattice_id_jump_pair.second;
    atom_id = lattice_id_jump_pair.first;
  }
 
 
  // Df at initial vacancy position
  auto Df_initial = GetDf(config, vacancy_id);

  // Df at final vacancy position, which is atom's initial position
  auto Df_final = GetDf(config, atom_id);

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




