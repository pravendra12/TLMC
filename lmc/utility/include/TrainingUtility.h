#ifndef LMC_CONSTANT_INCLUDE_TRAININGUTILITY_H_
#define LMC_CONSTANT_INCLUDE_TRAININGUTILITY_H_

#include <iostream>
#include <filesystem>
#include <nlohmann/json.hpp>

#include "Config.h"
#include "ConfigEncoding.h"
#include "SymmetrySpglib.h"
#include "BasisSet.h"
#include "ConfigEncoding.h"
#include "CorrelationVector.h"
#include "GetEquivalentClusters.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;
namespace fs = std::filesystem;

VectorXd GetKRAEncoding(
    const Config &config,
    const pair<size_t, size_t> &latticeIdPair,
    BasisSet &atomicBasis,
    const size_t &maxBondOrder,
    const size_t &maxClusterSize,
    const size_t &maxBondOrderOfCluster,
    const unordered_map<size_t, Eigen::RowVector3d> &referenceLatticeIdHashmap,
    const vector<pair<Eigen::Matrix3d, Eigen::Vector3d>> &symmetryOperations);

// Returns local cluster expansion around a lattice site
// To define local environment around a site which is
// Used for defining the local vacancy formation energy
VectorXd GetLocalSiteClusterVector(
    const Config &config,
    const size_t &latticeId,
    BasisSet &atomicBasis,
    const size_t maxBondOrder,
    const size_t maxClusterSize,
    const vector<pair<vector<vector<size_t>>, LatticeClusterType>> &encodedOrbits);

/*
nebOutput/
├─ 0_path1/
│  ├─ relaxation/
│  ├─ structures/
│  │  ├─ relaxed/
│  │  └─ unrelaxed/
│  │      ├─ POSCAR_unrelaxed_initial
│  │      └─ POSCAR_unrelaxed_final
│  └─ vacancyMigration/
├─ 0_path2/
│  ├─ relaxation/
│  ├─ structures/
│  │  ├─ relaxed/
│  │  └─ unrelaxed/
│  │      ├─ POSCAR_unrelaxed_initial
│  │      └─ POSCAR_unrelaxed_final
│  └─ vacancyMigration/
├─ 1_path1/
│  └─ ...
├─ 1_path2/
│  └─ ...
└─ ... (other *_path* directories)
*/
json ExtractTrainingDataForNEB(const string &pathToNEBOutput);

json ExtractLocalCEData(const string &pathToLocalCEOutput);

void IterateDirectoryToExtractData();

#endif // LMC_CONSTANT_INCLUDE_TRAININGUTILITY_H_
