#ifndef LMC_SYMCE_INCLUDE_SYMMETRICCE_H_
#define LMC_SYMCE_INCLUDE_SYMMETRICCE_H_

#include <map>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include "Config.h"
#include "Orbit.hpp"
#include "Structure.hpp"
#include "OrbitList.hpp"
#include "ClusterSpace.hpp"
#include "SymmetricCEInterface.h"
#include "LocalOrbitListGenerator.hpp"
#include "InterfaceUtility.h"
#include "SymmetrySpglib.h"
#include "BasisSet.h"
#include "AtomClusterType.hpp"


using namespace std;
using namespace Eigen;

class SymmetricCE
{
public:
  SymmetricCE(
      const Config &supercellConfig,
      const Config &primitiveConfig,
      const vector<string> &allowedElements,
      const vector<double> &clusterCutoffs,
      const double symprec = 1e-5);

  vector<vector<vector<int>>> GetLocalOrbitsForLatticeSite(
      const size_t &latticeId);

  vector<vector<vector<int>>> GetLocalOrbitsEncoding();

  vector<double> GetLocalClusterVector(
      const Config &config,
      const vector<size_t> &canonicalSortedLatticeIdVector,
      const vector<vector<vector<int>>> localOrbitEncoding) const;

  vector<double> GetClusterVectorForConfig(
    const Config &config);

  void TestSymmetricCE(
    Config &config, 
    const vector<double> &eciVector, 
    const pair<size_t, size_t> &latticeIdPair);

private:
  Config supercellConfig_;
  SymmetricCEInterface symCEInterface_;

  double fractionalPositionTolerance_;

  shared_ptr<OrbitList> primitiveOrbitList_;
  shared_ptr<Structure> supercellStructure_;
  ClusterSpace clusterSpace_;
};

#endif // LMC_SYMCE_INCLUDE_SYMMETRICCE_H_
