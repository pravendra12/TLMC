/**************************************************************************************************
 * Copyright (c) 2023-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 7/18/23 4:13 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 7/18/23 4:21 PM                                                           *
 **************************************************************************************************/

#ifndef LMC_CE_INCLUDE_SYMMETRY_H_
#define LMC_CE_INCLUDE_SYMMETRY_H_

#include "Config.h"

class Symmetry {

};



void FindPointGroupSymmetry(const Config &config, std::vector<size_t> cluster);
void FindPointGroupSymmetryMM2(const Config &config, std::vector<size_t> cluster);
std::pair<std::vector<Eigen::Matrix3i>, std::vector<Eigen::Vector3d>> getSymmetryOperations();
void equivalence_mapping(const Config &cfg, const std::vector<size_t> &cluster, double symprec, std::vector<std::vector<int>> &equivalence_map);
void find_equivalent_atoms_with_custom_point_group();
void FindPointGroupSymmetryCustom(const Config &config, std::vector<std::pair<size_t, Eigen::RowVector3d>> cluster);
std::pair<std::vector<Eigen::Matrix3d>, std::vector<Eigen::Vector3d>> GetRotationAndTranslationMatrix(const Config &config, const std::vector<size_t> &cluster);
void equivalent_custom(const Config &config);
void FindEquivalentClusters(const Config &config, const std::vector<std::vector<size_t>> &pair_clusters);
std::vector<int> SymmetricallySortLatticeIDs(const Config &config, const std::vector<size_t> &lattice_ids );


#endif // LMC_CE_INCLUDE_SYMMETRY_H_
