#include "Config.h"
#include "PotentialEnergyEstimator.h"
#include "VacancyMigrationPredictor.h"
#include "SymmetrySpglib.h"
#include "EnergyPredictor.h"
#include "TestingFunctions.h"
#include "ClusterExpansionParameters.h"
#include "BasisSet.h"
#include "AtomClusterType.hpp"
#include "CorrelationVector.h"

#include <filesystem>
#include <nlohmann/json.hpp>

namespace fs = std::filesystem;
using json = nlohmann::json;

int main()
{
  const vector<double> cutoffs = {3, 4, 5};

  auto cfgInitial = Config::ReadPoscar("/home/pravendra3/Documents/nebOutput/nebOutput/0_path2/structures/unrelaxed/POSCAR_unrelaxed_initial");
  auto cfgFinal = Config::ReadPoscar("/home/pravendra3/Documents/nebOutput/nebOutput/0_path2/structures/unrelaxed/POSCAR_unrelaxed_final");

  cfgInitial.UpdateNeighborList(cutoffs);
  cfgFinal.UpdateNeighborList(cutoffs);


  ClusterExpansionParameters ceParams("/home/pravendra3/Documents/LatticeMonteCarlo-eigen/script/coefficientFile_MoTa.json");

  EnergyPredictor energyPredictor(ceParams, cfgInitial);

  cout << "Energy Initial: " << energyPredictor.ComputeEnergyOfConfig(cfgInitial) << endl;
  cout << "Energy Final: " << energyPredictor.ComputeEnergyOfConfig(cfgFinal) << endl;

  cout << "dE: " << energyPredictor.ComputeEnergyOfConfig(cfgFinal) - energyPredictor.ComputeEnergyOfConfig(cfgInitial) << endl;


  auto vacancyId = cfgInitial.GetVacancyLatticeId();


  for (auto id : cfgInitial.GetNeighborLatticeIdVectorOfLattice(vacancyId, 1))
  {
    cout << energyPredictor.GetDeSwap(cfgInitial, make_pair(vacancyId, id)) << " : " << id << endl;


  }

  cout << cfgFinal.GetVacancyLatticeId() << endl;




}

/*
VectorXd GetKRAEncoding(

    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair,
    const size_t maxBondOrder,
    const unordered_map<size_t, RowVector3d> &canonicalReferenceMap,
    const vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>> &symmetryOperations,
    BasisSet &atomicBasis,
    const vector<pair<vector<vector<size_t>>,
                LatticeClusterType>> &equivalentClustersEncoding

)

{
  size_t id1 = latticeIdJumpPair.first;
  size_t id2 = latticeIdJumpPair.second;

  pair<size_t, size_t> canonicalJumpPair = (id1 < id2)
                                               ? std::make_pair(id1, id2)
                                               : std::make_pair(id2, id1);

  auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForPair(
      config,
      canonicalJumpPair,
      maxBondOrder,
      canonicalReferenceMap,
      symmetryOperations);

  VectorXd correlationVector = GetCorrelationVector(
      config,
      atomicBasis,
      canonicalSortedLatticeIds,
      equivalentClustersEncoding);

  return correlationVector;
}

VectorXd GetConfigEncoding(const Config &config,
                           size_t maxBondOrder,
                           BasisSet &atomicBasis,
                           const vector<pair<vector<vector<size_t>>,
                                       LatticeClusterType>> &equivalentClustersEncoding)
{
  // Assume all correlation vectors are the same size (depends only on basis + clusters)
  bool first = true;
  VectorXd totalCorrelation;

  for (size_t latticeId = 0; latticeId < config.GetNumLattices(); latticeId++)
  {
    // Get canonical cluster of sites involving this latticeId
    auto canonicalSortedLatticeIds = GetCanonicalSortedSitesForSite(
        config,
        latticeId,
        maxBondOrder);

    // Include the central site for consistency with encoding
    canonicalSortedLatticeIds.emplace_back(latticeId);

    // Compute correlation vector for this local cluster
    VectorXd correlationVector = GetCorrelationVector(
        config,
        atomicBasis,
        canonicalSortedLatticeIds,
        equivalentClustersEncoding);

    // Initialize accumulator on first pass
    if (first)
    {
      totalCorrelation = VectorXd::Zero(correlationVector.size());
      first = false;
    }

    // Add contribution from this site
    totalCorrelation += correlationVector;
  }

  return totalCorrelation;
}

int main()
{

  vector<double> cutoffs = {3, 4, 5};

  auto vacancy = Element("X");

  const string basisType = "Chebyshev";

  const size_t maxBondOrderCE = 3;
  const size_t maxClusterSizeCE = 3;

  const size_t maxBondOrderKRA = 2;
  const size_t maxClusterSizeKRA = 3;
  const size_t maxBondOrderOfClusterKRA = 3;

  set<Element> elementSetCE = {Element("Mo"), Element("Ta"), vacancy};
  set<Element> elementSetKRA = {Element("Mo"), Element("Ta")};

  BasisSet atomicBasisCE(
      elementSetCE,
      basisType,
      true);

  auto config = Config::GenerateSupercell(4, 3.2, "Mo", "BCC");
  config.UpdateNeighborList(cutoffs);

  bool isValidConfig = (config.GetNeighborLatticeIdVectorOfLattice(0, 1).size() == 8 &&
                        config.GetNeighborLatticeIdVectorOfLattice(0, 2).size() == 6 &&
                        config.GetNeighborLatticeIdVectorOfLattice(0, 3).size() == 12);

  cout << "isValidConfig: " << isValidConfig << endl;

  const auto eqClustersEncodingCE = GetEquivalentClustersEncoding(
      config,
      maxBondOrderCE,
      maxClusterSizeCE,
      true);

  const Vector3d jumpDirection(1, 1, 1);

  const unordered_map<size_t, RowVector3d> canonicalReferenceMap = GetCenteredNeighborsAlongJumpDirection(
      config,
      maxBondOrderKRA,
      jumpDirection);

  const auto eqClustersEncodingKRA = GetEquivalentClustersEncoding(
      config,
      maxBondOrderKRA,
      maxBondOrderOfClusterKRA,
      maxClusterSizeKRA,
      canonicalReferenceMap,
      true);

  BasisSet atomicBasisKRA(
      elementSetKRA,
      basisType,
      true);

  vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>> symmetryOperations = GetSymmetryOperations(config);

  // string pathNebOutput = "/home/pravendra3/Documents/nebOutput/nebOutput";

  fs::path pathNebOutput("/home/pravendra3/Documents/nebOutput/nebOutput/");

  json nebOutput;

  nebOutput["structureId"] = {};

  nebOutput["migratingElement"] = {};
  nebOutput["ceEncodingInitial"] = {};
  nebOutput["ceEncodingFinal"] = {};
  nebOutput["kraEncoding"] = {};

  for (const auto &entry : fs::directory_iterator(pathNebOutput))
  {
    cout << entry.path().filename() << endl;

    for (const auto &subEntry : fs::directory_iterator(entry))
    {

      if (subEntry.path().filename() == "structures")
      {
        cout << "\t" << subEntry.path().filename() << endl;

        string initialPoscarPath = subEntry.path().string() + "/unrelaxed/POSCAR_unrelaxed_initial";
        string finalPoscarPath = subEntry.path().string() + "/unrelaxed/POSCAR_unrelaxed_final";

        auto cfgInitial = Config::ReadPoscar(initialPoscarPath);
        auto cfgFinal = Config::ReadPoscar(finalPoscarPath);

        cfgInitial.UpdateNeighborList(cutoffs);
        cfgFinal.UpdateNeighborList(cutoffs);

        bool initialOk = (cfgInitial.GetNeighborLatticeIdVectorOfLattice(0, 1).size() == 8 &&
                          cfgInitial.GetNeighborLatticeIdVectorOfLattice(0, 2).size() == 6 &&
                          cfgInitial.GetNeighborLatticeIdVectorOfLattice(0, 3).size() == 12);

        bool finalOk = (cfgFinal.GetNeighborLatticeIdVectorOfLattice(0, 1).size() == 8 &&
                        cfgFinal.GetNeighborLatticeIdVectorOfLattice(0, 2).size() == 6 &&
                        cfgFinal.GetNeighborLatticeIdVectorOfLattice(0, 3).size() == 12);

        if (initialOk && finalOk)
        {
          // Cluster Expansion
          pair<size_t, size_t> jumpPair = {cfgInitial.GetVacancyLatticeId(), cfgFinal.GetVacancyLatticeId()};
          cout << "\t" << "Jump Pair: { " << jumpPair.first << ", " << jumpPair.second << " }" << endl;

          cout << "\t" << "Distance Order: " << cfgInitial.GetDistanceOrder(jumpPair.first, jumpPair.second) << endl;

          string migratingElement = cfgInitial.GetElementOfLattice(jumpPair.second).GetElementString();

          cout << "\t" << "Migrating Element: " << migratingElement << endl;

          VectorXd ceInitialEncoding = GetConfigEncoding(cfgInitial, maxBondOrderCE, atomicBasisCE, eqClustersEncodingCE);
          VectorXd ceFinalEncoding = GetConfigEncoding(cfgFinal, maxBondOrderCE, atomicBasisCE, eqClustersEncodingCE);

          VectorXd lceEncodingInitial = GetKRAEncoding(
              cfgInitial,
              jumpPair,
              maxBondOrderKRA,
              canonicalReferenceMap,
              symmetryOperations,
              atomicBasisKRA,
              eqClustersEncodingKRA);

          VectorXd lceEncodingFinal = GetKRAEncoding(
              cfgFinal,
              jumpPair,
              maxBondOrderKRA,
              canonicalReferenceMap,
              symmetryOperations,
              atomicBasisKRA,
              eqClustersEncodingKRA);

          if (lceEncodingFinal == lceEncodingInitial)
          {
            std::cout << "\tLCE are exactly equal" << endl;
          }

          string strId = entry.path().filename();

          nebOutput["structureId"].emplace_back(strId);
          nebOutput["migratingElement"].emplace_back(migratingElement);
          nebOutput["ceEncodingInitial"].emplace_back(ceInitialEncoding);
          nebOutput["ceEncodingFinal"].emplace_back(ceFinalEncoding);
          nebOutput["kraEncoding"].emplace_back(lceEncodingFinal);
        }
      }
    }
  }

  std::ofstream out("/home/pravendra3/Documents/nebOutput/nebOutputV2.json");
  out << nebOutput.dump(4); // pretty print with indent = 4
  out.close();
}
*/

/*
int main()
{path
  const vector<double> cutoffs = {3.3, 4.7, 5.6};
  auto cfg = Config::GenerateSupercell(5, 3.4, "W", "BCC");
  cfg.UpdateNeighborList(cutoffs);

  for (int i = 0; i < cfg.GetNumAtoms() / 3; i++)
  {
    cfg.SetElementOfLattice(i, Element("Ta"));
  }

  /*
  auto centralId = cfg.GetCentralAtomLatticeId();
  auto nnId = cfg.GetNeighborLatticeIdVectorOfLattice(centralId, 1)[0];

  Vector3d jumpDirection(1, 1, 1);
  auto referenceLatticeIdHashmap = GetCenteredNeighborsAlongJumpDirection(
      cfg, 2, jumpDirection);

  auto neighourLatticeIdVector = cfg.GetNeighborLatticeIdVectorOfLattice(0, 1);
  neighourLatticeIdVector.emplace_back(0);

  unordered_set<size_t> neighourLatticeIdSet(neighourLatticeIdVector.begin(),
                                             neighourLatticeIdVector.end());

  auto bccSymmetryOperations = GetSymmetryOperations(
      cfg, neighourLatticeIdSet, true);

  /*
  for (auto id : cfg.GetNeighborLatticeIdVectorOfLattice(137, 1))
  {

    cout << endl;
    cout << endl;
    cout << "Pair : 0-" << id << endl;
    cout << cfg.GetRelativeDistanceVectorLattice(137, id).transpose() << endl;

    auto sortedVectorOfPair = GetCanonicalSortedSitesForPair(
        cfg, make_pair(137, id), 2, referenceLatticeIdHashmap, bccSymmetryOperations);

    cout << "Sorted Vector: ";
    print1DVector(sortedVectorOfPair);
    cout << "Encoding: {";
    for (int i = 0; i < sortedVectorOfPair.size(); i++)
    {
      cout << i << ", ";
    }
    cout << " }" << endl;

    map<size_t, size_t> latticeIdToIndexMap;

    size_t idx = 0;
    for (auto id : sortedVectorOfPair)
    {
      latticeIdToIndexMap[id] = idx;
      idx++;
    }

    unordered_set<size_t> nnSet(sortedVectorOfPair.begin(), sortedVectorOfPair.end());

    // For dE

    vector<size_t> cluster = {137};

    auto neighbours = cfg.GetNeighborLatticeIdsUpToOrder(137, 1);
    neighbours.emplace_back(137);
    unordered_set<size_t> neighboursSet(neighbours.begin(), neighbours.end());

    auto latticeClusterCentered = FindAllLatticeClusters(cfg, 3, 1, cluster);

    auto siteCenteredEquivalentCluster = GetEquivalentClusters(
        cfg,
        neighboursSet,
        latticeClusterCentered, 1e-5, true);

    vector<size_t> cluster2 = {id};

    auto neighbours2 = cfg.GetNeighborLatticeIdsUpToOrder(id, 1);
    neighbours2.emplace_back(id);
    unordered_set<size_t> neighboursSet2(neighbours2.begin(), neighbours2.end());

    /*
    auto neighboursSet2 = cfg.GetNeighboringLatticeIdSetOfPair(make_pair(137, id), 2);
    neighboursSet2.insert(137);
    neighboursSet2.insert(id);


    auto latticeClusterCentered2 = FindAllLatticeClusters(cfg, 3, 1, cluster2);

    auto siteCenteredEquivalentCluster2 = GetEquivalentClusters(
        cfg,
        neighboursSet2,
        latticeClusterCentered2, 1e-5, true);
    /*
     cout << "Combined Cluster Based on Symmetry: " << endl;

     size_t countTriplets = 0;

     // Merge the two get the combined clusters
     std::vector<std::set<std::vector<size_t>>> combinedClusters;

     for (int i = 0; i < siteCenteredEquivalentCluster.size(); i++)
     {
       // set
       auto set1 = siteCenteredEquivalentCluster[i];
       auto set2 = siteCenteredEquivalentCluster2[i];

       set1.insert(set2.begin(), set2.end());

       combinedClusters.emplace_back(set1);
     }

     set<vector<size_t>> tripletSet1;

     for (size_t cid = 0; cid < combinedClusters.size(); ++cid)
     {
       cout << "Equivalence Class " << cid << ":\n";
       for (const auto &t : combinedClusters[cid])
       {
         cout << " Cluster: ";
         for (auto id : t)
           cout << id << " ";

         auto clusterType = IdentifyLatticeClusterType(cfg, t);
         cout << clusterType << "\n";

         if (t.size() == 3)
         {
           countTriplets++;
           tripletSet1.insert(t);
         }
       }
       cout << "\n";
     }

     cout << "Total number of triplets: " << countTriplets << endl;

     vector<size_t> test = {137, id};
     auto combinedClusterCe = FindAllLatticeClusters(
         cfg, 3, 1, test);

     cout << "All Combined Clusters : " << endl;

     std::vector<std::pair<LatticeClusterType, std::vector<size_t>>> sortedClusters;

     for (const auto &cluster : combinedClusterCe)
     {
       sortedClusters.emplace_back(
           cluster.GetClusterType(),
           cluster.GetLatticeIdVector());
     }
     // Sort by cluster type, then by lattice IDs
     std::sort(sortedClusters.begin(), sortedClusters.end(),
               [](const auto &a, const auto &b)
               {
                 if (a.first != b.first)
                   return a.first < b.first; // sort by type first
                 return a.second < b.second; // then lexicographically by IDs
               });

     // Print sorted clusters
     set<vector<size_t>> tripletSet2;

     countTriplets = 0;
     for (const auto &[type, latticeIds] : sortedClusters)
     {
       std::cout << "{ ";
       for (auto id : latticeIds)
       {
         std::cout << id << " ";
       }
       std::cout << "} " << type << "\n";

       if (latticeIds.size() == 3)
       {
         countTriplets += 1;
         tripletSet2.insert(latticeIds);
       }
     }

     cout << "Size of tripletSet1: " << tripletSet1.size() << endl;
     cout << "Size of tripletSet2: " << tripletSet2.size() << endl;

     cout << "Total number of triplets: " << countTriplets << endl;

     std::set<std::vector<size_t>> symDiff;

     std::set_intersection(
         siteCenteredEquivalentCluster2[1].begin(), siteCenteredEquivalentCluster2[1].end(),
         siteCenteredEquivalentCluster2[3].begin(), siteCenteredEquivalentCluster2[3].end(),
         std::inserter(symDiff, symDiff.begin()));

     cout << "Symmetric difference (elements NOT in intersection):\n";

     for (const auto &v : symDiff)
     {
       for (auto id : v)
         std::cout << id << " ";
       std::cout << "\n";
     }
     */

//

// In KRA no need to include the original sites
/*
auto latticeClusterKRA = FindClustersWithinAllowedSites(
    cfg, 3, 3, sortedVectorOfPair);

auto kraEquivalentCluster = GetEquivalentClusters(
    cfg,
    nnSet,
    latticeClusterKRA, 1e-5, true);

cout << "Equivalent Cluster Encoding Site Centered: " << endl;

vector<vector<vector<size_t>>> equivalentClusterGlobal;
for (size_t cid = 0; cid < kraEquivalentCluster.size(); ++cid)
{
  cout << "Equivalence Class " << cid << ":\n";
  vector<vector<size_t>> equivalentClusters;
  for (const auto &t : kraEquivalentCluster[cid])
  {
    vector<size_t> encodedVector;
    for (auto id : t)
      encodedVector.emplace_back(latticeIdToIndexMap.at(id));
    sort(encodedVector.begin(), encodedVector.end());
    equivalentClusters.emplace_back(encodedVector);
  }
  sort(equivalentClusters.begin(), equivalentClusters.end());
  print2DVector(equivalentClusters);

  equivalentClusterGlobal.emplace_back(equivalentClusters);
  cout << "\n";
}


sort(equivalentClusterGlobal.begin(), equivalentClusterGlobal.end());
print3DVector(equivalentClusterGlobal);

cout << id << endl;
*/

/*
size_t id = 137;

vector<size_t> cluster2 = {id};

auto neighbours2 = cfg.GetNeighborLatticeIdsUpToOrder(id, 1);
unordered_set<size_t> neighboursSet2(neighbours2.begin(), neighbours2.end());

/*
auto neighboursSet2 = cfg.GetNeighboringLatticeIdSetOfPair(make_pair(137, id), 2);
neighboursSet2.insert(137);
neighboursSet2.insert(id);


auto latticeClusterCentered2 = FindAllLatticeClusters(cfg, 3, 1, cluster2);

auto nnVector = cfg.GetNeighborLatticeIdsUpToOrder(0, 1);
unordered_set<size_t> nnSet(nnVector.begin(), nnVector.end());

auto symmetryOperations = GetSymmetryOperations(cfg, nnSet, true);

neighboursSet2.insert(id);

GetEquivalentClusters(
    cfg,
    neighboursSet2,
    latticeClusterCentered2,
    1e-5,
    true);

auto nnSetPair = cfg.GetNeighboringLatticeIdSetOfPair(make_pair(137, 112), 2);

auto symOps = GetSymmetryOperations(cfg, nnSetPair, true);

auto rfMap = GetCenteredNeighborsAlongJumpDirection(cfg, 2, jumpDirection);

auto clustersEncoding = GetEquivalentClustersEncoding(
    cfg, 2, 2, rfMap, true);

for (auto id : cfg.GetNeighborLatticeIdVectorOfLattice(0, 1))
{
  std::cout << "0-" << id << std::endl
            << std::endl;

  // Get canonical sorted lattice sites for this pair
  auto sortedLatticeIds = GetCanonicalSortedSitesForPair(cfg, std::make_pair(0, id), 2, rfMap, bccSymmetryOperations);

  // Encode the lattice IDs
  std::map<size_t, size_t> latticeIdHashMap; // or use unordered_map if faster
  std::vector<size_t> encodedCluster;
  int i = 0;
  for (auto latticeId : sortedLatticeIds)
  {
    latticeIdHashMap[latticeId] = i;
    i++;
  }

  unordered_set<size_t> latticeIdSet(sortedLatticeIds.begin(), sortedLatticeIds.end());

  GetSymmetryOperations(cfg, latticeIdSet, true, latticeIdHashMap);

  // break;
}

// KRA Encoding
auto sortedLatticeIds = GetCanonicalSortedSitesForPair(cfg, std::make_pair(0, 25), 2, rfMap, bccSymmetryOperations);

auto atomVector = cfg.GetAtomVector();

set<Element> elementSet(atomVector.begin(), atomVector.end());
string basisType = "Occupation";

auto encoding = GetLocalEnvironmentEncoding(cfg, elementSet, basisType, sortedLatticeIds, clustersEncoding);

cout << encoding << endl;

// for computing dE

print1DVector(cfg.GetNeighborLatticeIdsUpToOrder(0, 3));
print1DVector(cfg.GetNeighborLatticeIdVectorOfLattice(0, 1));
print1DVector(cfg.GetNeighborLatticeIdVectorOfLattice(0, 2));
print1DVector(cfg.GetNeighborLatticeIdVectorOfLattice(0, 3));

GetEquivalentClustersEncoding(
    cfg,
    2,
    2,
    true);

*/
/*
set<Element> elementSet = {Element("Ta"), Element("W"), Element("X")};
string basisType = "Occupation";

for (size_t i : {123, 167})
{
  cout << "----------------------------" << endl;
  cout << "\nSite " << i << endl;

  auto nnLatticeIds = GetCanonicalSortedSitesForSite(cfg, i, 3);

  nnLatticeIds.emplace_back(i);
  unordered_set<size_t> nnLatticeIdSet(nnLatticeIds.begin(), nnLatticeIds.end());

  auto equivalentClustersEncoding = GetEquivalentClustersEncoding(
      cfg,
      3,
      3,
      false);
}

cout << cfg.GetElementOfLattice(123) << endl;
cout << cfg.GetElementOfLattice(167) << endl;

vector<string> elementStrings = ReadStringParametersFromJson("/home/pravendra3/Documents/LatticeMonteCarlo-eigen/bin/predictorFileKRATesting.json",
                                                             "ce", "elementSet");

print1DVector(elementStrings);


TestEnergyChangePredictor(cfg, predictorFilename);

ClusterExpansionParameters ceParams(predictorFilename, true);

LruCache<size_t, VectorXd> lruTest(20);
VectorXd someVector(4); // size 4
someVector << 0, 1, 2, 3;

Element e("Ta");
size_t elementKey = boost::hash<Element>()(e);

lruTest.Add(elementKey, someVector);

cout << elementKey << endl;

size_t elementKey2 = boost::hash<Element>()(e);

cout << elementKey2 << endl;

cout << lruTest.GetSizeOfCacheList() << endl;

AtomClusterType atomClusterType(Element("Ta"), Element("Al"));

auto atomClusterHashKey = boost::hash<AtomClusterType>()(atomClusterType);

cout << atomClusterHashKey << endl;

AtomClusterType atomClusterType2;

auto atomClusterHashKey2 = boost::hash<AtomClusterType>()(atomClusterType2);

cout << atomClusterHashKey2 << endl;

*/
/*
  auto eqClustersEncoding = GetEquivalentClustersEncoding(
      cfg,
      3,
      3,
      true);

  auto canSortedSites = GetCanonicalSortedSitesForSite(
      cfg,
      0,
      3);

  print1DVector(canSortedSites);
  canSortedSites.emplace_back(0);

  print1DVector(canSortedSites);


  auto atomVector = cfg.GetAtomVector();
  set<Element> elementSetOfCfg(atomVector.begin(), atomVector.end());

  BasisSet atomicBasis(
      elementSetOfCfg,
      "Chebyshev",
      true);

  auto featureVector = GetCorrelationVector(
      cfg,
      atomicBasis,
      canSortedSites,
      eqClustersEncoding);

  cout << featureVector.transpose() << endl;

  string predictorFilename = "/home/pravendra3/Documents/LatticeMonteCarlo-eigen/bin/predictorFileKRATesting.json";

  ClusterExpansionParameters ceParams(predictorFilename);

  EnergyPredictor energyPredictor(ceParams, cfg);

  energyPredictor.ComputeEnergyOfSite(cfg, 100);

  // energyPredictor.ComputeEnergyOfSite(cfg, 0);

  energyPredictor.GetDeMigration(cfg, make_pair(0, cfg.GetNeighborLatticeIdVectorOfLattice(0, 1)[0]));
}
*/
