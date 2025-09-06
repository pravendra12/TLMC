#ifndef LMC_SYMCE_INCLUDE_SYMMETRICCEINTERFACE_H_
#define LMC_SYMCE_INCLUDE_SYMMETRICCEINTERFACE_H_

#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include <string>
#include "Config.h"
#include <filesystem>
#include "InterfaceUtility.h"
#include "Structure.hpp"
#include "LatticeSite.hpp"

using namespace std;
namespace py = pybind11;
namespace fs = std::filesystem;

class SymmetricCEInterface
{
public:
  // config : Generally a small supercell or a primitive structure
  // allowedElements : Elements which are allowed to occupy sites of config
  // clusterCuttoffs : Cutoffs for {pair, triplet, quadruplet} clusters
  SymmetricCEInterface(
      const Config &config,
      const vector<string> &allowedElements,
      const vector<double> &clusterCutoffs,
      const double symprec = 1e-5);

  // Inputs Required to initialize OrbitList and ClusterSpace

  // Returns primitve structure
  Structure GetPrimitiveStructure();

  // Returns matrix of equivalent sites
  vector<vector<LatticeSite>> GetMatrixOfEquivalentSites();

  // Returns neighbour lists
  vector<vector<vector<LatticeSite>>> GetNeighbourLists();

  // Returns tolerance
  double GetPositionTolerance();

  // Returns fractional position tolerance
  double GetFractionalPositionTolerance();

private:
  // Convert Config to ASE Atoms
  py::object ConvertConfigToASEAtoms(
      py::module_ &aseAtoms,
      const Config &config);

  py::scoped_interpreter pyGuard_;

  // Utility Module
  py::module_ pyStandalone_;

  // ASE Atoms Object
  py::object pyAtomsObject_;

  // Cluster Cutoffs
  py::list pyCutoffs_;

  // Allowed chemical species
  py::list pyAllowedSpecies_;

  // Allowed Atomic Numbers // vector<vector<int>>
  vector<vector<int>> allowedAtomicNumbersVector_;

  // Tolerance
  double symprec_;
};

#endif // LMC_SYMCE_INCLUDE_SYMMETRICCEINTERFACE_H_
