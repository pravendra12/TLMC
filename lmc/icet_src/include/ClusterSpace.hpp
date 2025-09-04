#pragma once

#define _USE_MATH_DEFINES
#include <cmath>

#include "LocalOrbitListGenerator.hpp"
#include "OrbitList.hpp"
#include "Structure.hpp"
#include "VectorOperations.hpp"

/**
@brief This class handles the cluster space.
@details It provides functionality for setting up a cluster space, calculating
cluster vectors as well as retrieving various types of associated information.
*/
class ClusterSpace
{
public:
  /// Constructor.
  ClusterSpace(){};
  ClusterSpace(std::shared_ptr<OrbitList>, const double, const double);

  /// Returns the entire primitive orbit list.
  const OrbitList &getPrimitiveOrbitList() const { return *_primitiveOrbitList; }

  /// Returns the primitive structure.
  const Structure &primitiveStructure() const { return _primitiveOrbitList->structure(); }

  /// Returns the cluster space size, i.e. the length of a cluster vector.
  size_t size() const;

  /// Returns the mapping between atomic numbers and the internal species enumeration scheme for each site.
  const std::vector<std::unordered_map<int, int>> &getSpeciesMaps() const { return _speciesMaps; }

  /// Returns the cluster vector corresponding to the input structure.
  std::vector<double> getClusterVector(const Structure &, const double) const;

  /// Returns the cluster vector given the orbit list and a structure.
  const std::vector<double> getClusterVectorFromOrbitList(
    const OrbitList &, const std::shared_ptr<Structure>, const int flipIndex = -1, const int newOccupation = -1) const;

  const std::vector<double> getClusterVectorFromOrbitListCheck(
    const OrbitList &orbitList, const std::shared_ptr<Structure> supercell, const int flipIndex = -1, const int newOccupation = -1) const;

  /// Returns the cluster product.
  double evaluateClusterProduct(const std::vector<int> &, const std::vector<int> &, const std::vector<int> &, const std::vector<int> &, const std::vector<int> &) const;

  /// Returns the default cluster function.
  double evaluateClusterFunction(const int, const int, const int) const;

  /// Merge orbits.
  void mergeOrbits(const int index1, const int index2) { _primitiveOrbitList->mergeOrbits(index1, index2); }

private:
  /// Primitive orbit list based on the structure and the cutoffs.
  std::shared_ptr<OrbitList> _primitiveOrbitList;

  /// Map between atomic numbers and the internal species enumeration scheme for each site in the primitive structure.
  std::vector<std::unordered_map<int, int>> _speciesMaps;
};
