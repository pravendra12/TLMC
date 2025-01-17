
#ifndef LMC_CE_INCLUDE_SYMMETRYCUSTOM_H_
#define LMC_CE_INCLUDE_SYMMETRYCUSTOM_H_

#include "Config.h"
#include "eigen3/Eigen/Dense"


std::vector<size_t> GetSymmetricallySortedLatticeVectorMMM(const Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair, const size_t max_bond_order);

std::vector<size_t> GetSymmetricallySortedLatticeVectorMM2(const Config &config, const std::pair<size_t, size_t> &lattice_id_jump_pair, const size_t max_bond_order);

#endif //LMC_CE_INCLUDE_SYMMETRYCUSTOM_H_

