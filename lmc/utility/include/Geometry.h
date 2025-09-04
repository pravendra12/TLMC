#ifndef _LMC_UTILITY_INCLUDE_GEOMETRY_H_
#define _LMC_UTILITY_INCLUDE_GEOMETRY_H_

#include <Eigen/Dense>
#include <vector>
#include <string>
#include <spglib.h>
#include "Structure.hpp"
#include "LatticeSite.hpp"
#include "Element.hpp"
#include <set>
#include <map>


using namespace std;

Eigen::MatrixXd get_scaled_positions(const Eigen::MatrixXd &positions, const Eigen::Matrix3d &cell, bool wrap, const std::array<bool, 3> &pbc);

Eigen::MatrixXd fractional_to_cartesian(const Structure & structure, const Eigen::MatrixXd & frac_positions);
std::vector<std::string> get_wyckoff_sites(const Structure &structure, const std::vector<std::vector<std::string>> &map_occupations, double symprec, bool include_representative_atom_index);

Structure get_primitive_structure(const Structure &structure, bool no_idealize, bool to_primitive, double symprec);
std::pair<Structure, std::vector<std::vector<std::string>>> get_occupied_primitive_structure(const Structure &structure, const std::vector<std::vector<std::string>> &allowed_species, double symprec);
#endif // _LMC_UTILITY_INCLUDE_GEOMETRY_H_
