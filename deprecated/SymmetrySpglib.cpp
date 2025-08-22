#include "SymmetrySpglib.h"
#include "UnionFind.h"
#include "ClusterExpansion.h"
#include "PrintUtility.h"

void get_bcc_symmetry_operations(
    Config &cfg,
    vector<size_t> &lce,
    unordered_set<Matrix3d, Matrix3dHash> &uniqueRotations,
    unordered_set<Vector3d, Vector3dHash2> &uniqueTranslations,
    vector<pair<Matrix3d, Vector3d>> &symmetryOperations)
{
  const double symprec = 1e-5;

  // BCC primitive cell
  double lattice[3][3] = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}};

  int numAtoms = lce.size();

  double (*positions)[3] = new double[numAtoms][3];
  int *types = new int[numAtoms];

  for (int i = 0; i < numAtoms; i++)
  {
    VectorXd relativePosition = cfg.GetRelativePositionOfLattice(lce[i]);

    positions[i][0] = relativePosition(0);
    positions[i][1] = relativePosition(1);
    positions[i][2] = relativePosition(2);

    types[i] = 1;
  }

  SpglibDataset *dataset = spg_get_dataset(lattice, positions, types, numAtoms, symprec);

  if (!dataset)
  {
    cerr << "Failed to get spglib dataset." << endl;
    return;
  }

  cout << "Space group: " << dataset->international_symbol
       << " (#" << dataset->spacegroup_number << ")" << endl;

  cout << "Point group: " << dataset->pointgroup_symbol << endl;

  for (int i = 0; i < dataset->n_operations; ++i)
  {
    Matrix3d rot;
    Vector3d trans;

    for (int r = 0; r < 3; ++r)
    {
      trans(r) = dataset->translations[i][r];
      for (int c = 0; c < 3; ++c)
      {
        rot(r, c) = dataset->rotations[i][r][c];
      }
    }

    uniqueRotations.insert(rot);
    uniqueTranslations.insert(trans);
    pair<Matrix3d, Vector3d> symOpPair = {rot, trans};
    symmetryOperations.emplace_back(symOpPair);
  }

  // cout << "\nUnique Rotation Matrices: " << uniqueRotations.size() << endl;
  // for (const auto &rot : uniqueRotations)
  // {
  //   cout << rot << "\n---\n";
  // }

  // cout << "\nUnique Translation Vectors: " << uniqueTranslations.size() << endl;
  // for (const auto &trans : uniqueTranslations)
  // {
  //   cout << trans.transpose() << endl;
  // }

  // cout << "\nTransformation matrix to conventional cell:\n";
  // for (int i = 0; i < 3; ++i)
  // {
  //   // cout << "[ ";
  //   for (int j = 0; j < 3; ++j)
  //   {
  //     cout << dataset->transformation_matrix[i][j] << " ";
  //   }
  //   cout << "]" << endl;
  // }

  map<size_t, vector<size_t>> equivalence_map;

  for (int i = 0; i < dataset->n_atoms; i++)
  {
    size_t key = dataset->equivalent_atoms[i];
    size_t value = lce[i];

    if (equivalence_map.find(key) != equivalence_map.end())
    {
      equivalence_map[key].push_back(value); // use push_back instead of insert
    }
    else
    {
      equivalence_map[key] = {value};
    }
    // cout << value << " - " << key << endl;
  }

  for (auto group : equivalence_map)
  {
    cout << group.first << " : { ";

    for (auto id : group.second)
    {
      cout << id << " ";
    }
    cout << " } " << endl;
  }

  spg_free_dataset(dataset);
  delete[] positions;
}

struct Vector3dEqual
{
  bool operator()(const Eigen::Vector3d &a, const Eigen::Vector3d &b) const
  {
    return (a - b).norm() < 1e-6;
  }
};

#include <set>

// Helper to check if two clusters (vectors of lattice IDs) are equal up to ordering
bool clustersAreSame(const vector<size_t> &a, const vector<size_t> &b)
{
  if (a.size() != b.size())
    return false;
  auto a_sorted = a;
  std::sort(a_sorted.begin(), a_sorted.end());
  auto b_sorted = b;
  std::sort(b_sorted.begin(), b_sorted.end());
  return a_sorted == b_sorted;
}

// Function to apply symmetry operation (rot + trans) to cluster positions
vector<Vector3d> applySymOpToCluster(const vector<Vector3d> &positions,
                                     const Matrix3d &rot,
                                     const Vector3d &trans)
{
  vector<Vector3d> transformed;
  for (const auto &pos : positions)
  {
    Vector3d p = rot * pos + trans;
    // wrap fractional coordinates to [0,1)
    for (int i = 0; i < 3; i++)
    {
      p[i] = fmod(p[i], 1.0);
      if (p[i] < 0)
        p[i] += 1.0;
    }
    transformed.push_back(p);
  }
  return transformed;
}

// Function to find lattice IDs for cluster positions with tolerance
// Using your latticeIdToPositionMap and tolerance
vector<size_t> getLatticeIdsFromPositions(const vector<Vector3d> &positions,
                                          const unordered_map<size_t, Vector3d> &latticeIdToPos,
                                          double tol)
{
  vector<size_t> ids;
  for (const auto &pos : positions)
  {
    bool found = false;
    for (const auto &pair : latticeIdToPos)
    {
      if ((pos - pair.second).norm() < tol)
      {
        ids.push_back(pair.first);
        found = true;
        break;
      }
    }
    if (!found)
      throw std::runtime_error("Position not found within tolerance");
  }
  return ids;
}

void equivalent_clusters(Config &cfg, vector<size_t> lce)
{
  unordered_set<Matrix3d, Matrix3dHash> uniqueRotations;
  unordered_set<Vector3d, Vector3dHash2> uniqueTranslations;
  vector<pair<Matrix3d, Vector3d>> symOperations;

  get_bcc_symmetry_operations(cfg, lce, uniqueRotations, uniqueTranslations, symOperations);

  // Convert unordered_sets to vectors for easier iteration
  vector<Matrix3d> rotations(uniqueRotations.begin(), uniqueRotations.end());
  vector<Vector3d> translations(uniqueTranslations.begin(), uniqueTranslations.end());

  std::unordered_map<size_t, Vector3d> latticeIdToPositionMap;

  for (auto latticeId : lce)
  {
    latticeIdToPositionMap[latticeId] = cfg.GetRelativePositionOfLattice(latticeId);
  }

  auto centralId = cfg.GetCentralAtomLatticeId();
  vector<size_t> latticeCluster = {centralId};

  auto clusterAll = FindAllLatticeClusters(cfg, 2, 3, latticeCluster);

  vector<vector<size_t>> pairClusters;

  for (auto cluster : clusterAll)
  {
    auto latticeIds = cluster.GetLatticeIdVector();

    if (latticeIds.size() == 2)
    {
      pairClusters.emplace_back(latticeIds);
    }
  }

  const double tolerance = 1e-5;

  // Now group pairClusters into equivalence classes
  vector<set<size_t>> equivalenceClasses; // each set contains indices of pairClusters equivalent to each other
  vector<bool> visited(pairClusters.size(), false);

  for (size_t i = 0; i < pairClusters.size(); ++i)
  {
    if (visited[i])
      continue;

    set<size_t> eqClass;
    eqClass.insert(i);
    visited[i] = true;

    // Get positions of cluster i
    vector<Vector3d> clusterPos;
    for (auto id : pairClusters[i])
      clusterPos.push_back(latticeIdToPositionMap[id]);

    // Now compare with all other clusters j > i
    for (size_t j = i + 1; j < pairClusters.size(); ++j)
    {
      if (visited[j])
        continue;

      // Get positions of cluster j
      vector<Vector3d> clusterPosJ;
      for (auto id : pairClusters[j])
        clusterPosJ.push_back(latticeIdToPositionMap[id]);

      bool equivalent = false;

      // Check if j is symmetry equivalent to i by trying all symmetry operations
      for (size_t k = 0; k < rotations.size(); ++k)
      {
        try
        {
          vector<Vector3d> transformed = applySymOpToCluster(clusterPos, rotations[k], translations[k]);
          vector<size_t> transformedIds = getLatticeIdsFromPositions(transformed, latticeIdToPositionMap, tolerance);

          // Sort to canonicalize
          std::sort(transformedIds.begin(), transformedIds.end());
          auto jSorted = pairClusters[j];
          std::sort(jSorted.begin(), jSorted.end());

          if (transformedIds == jSorted)
          {
            equivalent = true;
            break;
          }
        }
        catch (std::runtime_error &)
        {
          // Position not found, try next symmetry op
        }
      }

      if (equivalent)
      {
        eqClass.insert(j);
        visited[j] = true;
      }
    }

    equivalenceClasses.push_back(eqClass);
  }

  // Print equivalence classes (cluster indices and their lattice IDs)
  int classId = 0;
  for (const auto &eqClass : equivalenceClasses)
  {
    std::cout << "Equivalence Class " << classId++ << ":\n";
    for (auto idx : eqClass)
    {
      std::cout << " Cluster " << idx << ": ";

      for (auto id : pairClusters[idx])
        std::cout << id << " ";
      std::cout << " " << cfg.GetDistanceOrder(pairClusters[idx][0], pairClusters[idx][1]) << "\n";
    }
    std::cout << "\n";
  }
}
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <cmath>
#include <algorithm>
#include <utility> // for std::pair

using namespace std;

#include <iostream>
#include <iomanip>

// Assuming Vector3d and Matrix3d are Eigen types or similar
// and apply_symmetry_op() and latticeIdToPositionMap, symOperations are defined

void debug_compare_these_two(
    const unordered_map<size_t, Vector3d> &latticeIdToPositionMap,
    const vector<pair<Matrix3d, Vector3d>> &symOperations,
    double tolerance = 1e-5)
{
  vector<size_t> clusterA = {112, 137, 118};
  vector<size_t> clusterB = {162, 137, 168};

  vector<Vector3d> posA, posB;
  for (auto id : clusterA)
    posA.push_back(latticeIdToPositionMap.at(id));
  for (auto id : clusterB)
    posB.push_back(latticeIdToPositionMap.at(id));

  std::cout << "Comparing clusters:\n A: ";
  for (auto id : clusterA)
    std::cout << id << " ";
  std::cout << "\n B: ";
  for (auto id : clusterB)
    std::cout << id << " ";
  std::cout << "\n\n";

  int opCount = 0;
  for (const auto &[rot, trans] : symOperations)
  {
    std::cout << "Symmetry operation #" << opCount++ << ":\n";

    std::cout << "Rotation matrix:\n";
    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
        std::cout << std::setw(10) << std::setprecision(5) << rot(i, j) << " ";
      std::cout << "\n";
    }

    std::cout << "Translation vector: ";
    for (int i = 0; i < 3; ++i)
      std::cout << std::fixed << std::setprecision(5) << trans[i] << " ";
    std::cout << "\n";

    // Apply symmetry op to clusterA positions
    vector<Vector3d> transformed;
    for (const auto &pos : posA)
    {
      Vector3d p = rot * pos + trans;
      for (int i = 0; i < 3; ++i)
      {
        p[i] = std::fmod(p[i], 1.0);
        if (p[i] < 0)
          p[i] += 1.0;
      }
      transformed.push_back(p);
    }

    std::cout << "Transformed positions of cluster A:\n";
    for (size_t i = 0; i < transformed.size(); ++i)
    {
      std::cout << "  Site " << clusterA[i] << ": ("
                << std::fixed << std::setprecision(5)
                << transformed[i][0] << ", "
                << transformed[i][1] << ", "
                << transformed[i][2] << ")\n";
    }

    // Check matches with cluster B (unordered)
    vector<bool> matched(posB.size(), false);
    bool allMatched = true;
    for (size_t i = 0; i < transformed.size(); ++i)
    {
      bool foundMatch = false;
      for (size_t j = 0; j < posB.size(); ++j)
      {
        if (!matched[j] && (transformed[i] - posB[j]).norm() < tolerance)
        {
          matched[j] = true;
          foundMatch = true;
          std::cout << "  Transformed site " << clusterA[i] << " matches cluster B site " << clusterB[j] << "\n";
          break;
        }
      }
      if (!foundMatch)
      {
        allMatched = false;
        std::cout << "  Transformed site " << clusterA[i] << " has NO match in cluster B\n";
      }
    }

    if (allMatched)
      std::cout << "=> Clusters ARE equivalent under this symmetry operation!\n";
    else
      std::cout << "=> Clusters NOT equivalent under this symmetry operation.\n";

    std::cout << "-------------------------------------------\n";
  }
}

// Helper: wrap fractional coordinate in [0,1)
inline double mod1(double x)
{
  x = std::fmod(x, 1.0);
  if (x < 0)
    x += 1.0;
  return x;
}

// Helper: apply symmetry op (rotation + translation) and wrap positions
vector<Vector3d> apply_symmetry_op(const vector<Vector3d> &positions,
                                   const Matrix3d &rot, const Vector3d &trans)
{
  vector<Vector3d> transformed;
  for (const auto &pos : positions)
  {
    Vector3d p = rot * pos + trans;
    for (int i = 0; i < 3; ++i)
      p[i] = mod1(p[i]);
    transformed.push_back(p);
  }
  return transformed;
}

vector<Vector3d> canonicalize_cluster(const vector<Vector3d> &positions)
{
  vector<Vector3d> canonical = positions; // copy

  // Sort by lex order: first x, then y, then z
  std::sort(canonical.begin(), canonical.end(), [](const Vector3d &a, const Vector3d &b)
            {
        if (a[0] != b[0]) return a[0] < b[0];
        if (a[1] != b[1]) return a[1] < b[1];
        return a[2] < b[2]; });

  return canonical;
}
/*
void equivalent_clusters_triplets(Config &cfg, vector<size_t> lce)
{
  double tolerance = 1e-5;
  unordered_set<Matrix3d, Matrix3dHash> uniqueRotations;
  unordered_set<Vector3d, Vector3dHash2> uniqueTranslations;
  vector<pair<Matrix3d, Vector3d>> symOperations;

  // Get symmetry operations for your BCC lattice
  get_bcc_symmetry_operations(cfg, lce, uniqueRotations, uniqueTranslations, symOperations);

  unordered_map<size_t, Vector3d> latticeIdToPositionMap;
  for (auto latticeId : lce)
    latticeIdToPositionMap[latticeId] = cfg.GetRelativePositionOfLattice(latticeId);

  auto centralId = cfg.GetCentralAtomLatticeId();
  vector<size_t> latticeCluster = {centralId};

  // Add centralId if missing
  if (find(lce.begin(), lce.end(), centralId) == lce.end())
    lce.emplace_back(centralId);

  // Find all clusters of size 3 within allowed sites
  auto clusterAll = FindClustersWithinAllowedSites(cfg, 3, 3, lce);

  // Collect only triplets containing centralId
  vector<vector<size_t>> tripletClusters;
  for (const auto &cluster : clusterAll)
  {
    auto latticeIds = cluster.GetLatticeIdVector();
    if (latticeIds.size() == 3 &&
        find(latticeIds.begin(), latticeIds.end(), centralId) != latticeIds.end())
    {
      tripletClusters.push_back(latticeIds);
    }
  }

  vector<bool> visited(tripletClusters.size(), false);
  vector<set<size_t>> equivalenceClasses;

  for (size_t i = 0; i < tripletClusters.size(); ++i)
  {
    if (visited[i])
      continue;

    set<size_t> eqClass;
    eqClass.insert(i);
    visited[i] = true;

    // Positions of cluster i
    vector<Vector3d> clusterPosI;
    for (auto id : tripletClusters[i])
      clusterPosI.push_back(latticeIdToPositionMap[id]);

    for (size_t j = i + 1; j < tripletClusters.size(); ++j)
    {
      if (visited[j])
        continue;

      vector<Vector3d> clusterPosJ;
      for (auto id : tripletClusters[j])
        clusterPosJ.push_back(latticeIdToPositionMap[id]);

      bool equivalent = false;

      // Check all symmetry operations
      for (const auto &[rot, trans] : symOperations)
      {
        auto canonicalCluster = canonicalize_cluster(clusterPosI);
        auto transformed = apply_symmetry_op(canonicalCluster, rot, trans);

        // Try to match transformed positions to lattice IDs
        vector<size_t> transformedIds;
        bool allFound = true;

        for (const auto &pos : transformed)
        {
          bool found = false;
          for (const auto &pair : latticeIdToPositionMap)
          {
            if ((pos - pair.second).norm() < tolerance)
            {
              transformedIds.push_back(pair.first);
              found = true;
              break;
            }
          }
          if (!found)
          {
            allFound = false;
            break;
          }
        }

        if (!allFound)
          continue;

        // Sort for comparison
        std::sort(transformedIds.begin(), transformedIds.end());
        auto jSorted = tripletClusters[j];
        std::sort(jSorted.begin(), jSorted.end());

        if (transformedIds == jSorted)
        {
          equivalent = true;
          break;
        }
      }

      if (equivalent)
      {
        eqClass.insert(j);
        visited[j] = true;
      }
      else
      {
        // Optional debug print
        // cout << "Clusters " << i << " and " << j << " NOT equivalent\n";
      }
    }

    equivalenceClasses.push_back(eqClass);
  }

  // Print equivalence classes
  int classId = 0;
  for (const auto &eqClass : equivalenceClasses)
  {
    cout << "Equivalence Class " << classId++ << ":\n";
    for (auto idx : eqClass)
    {
      cout << " Cluster " << idx << ": ";
      for (auto id : tripletClusters[idx])
        cout << id << " ";
      cout << " Distances: ";
      for (auto id : tripletClusters[idx])
        cout << cfg.GetDistanceOrder(centralId, id) << " ";
      cout << "\n";
    }
    cout << "\n";
  }

  // debug_compare_these_two(latticeIdToPositionMap, symOperations);


}
*/
#include <iostream>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <set>

using namespace std;
#include "UnionFind.h"

#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <queue>

using namespace std;

// NOTE: adapt types and helper functions (Matrix3d, Vector3d etc.) to your project.
static vector<size_t> canonicalize_triplet(vector<size_t> t)
{
  sort(t.begin(), t.end());
  return t;
}

static string triplet_to_string(const vector<size_t> &t)
{
  // small string key useful for printing or unordered_map keys if needed
  ostringstream oss;
  for (size_t i = 0; i < t.size(); ++i)
  {
    if (i)
      oss << "_";
    oss << t[i];
  }
  return oss.str();
}

// Quantize fractional coords to a stable string key (tunable grid)
static string frac_key(const Vector3d &v, double grid = 1e6)
{
  auto q = [&](double x) -> long long
  {
    double xx = x - floor(x); // bring to [0,1)
    // round to nearest grid cell
    return llround(xx * grid);
  };
  ostringstream oss;
  oss << q(v[0]) << "_" << q(v[1]) << "_" << q(v[2]);
  return oss.str();
}

// ---------- form_equivalence_classes (BFS) ---------------------------------
//
// Input: map<canonical_triplet, set<canonical_triplet>> equivalent_map
// Output: vector<set<canonical_triplet>> equivalence classes
// and mapping triplet -> classId
//
static vector<set<vector<size_t>>>
form_equivalence_classes(const map<vector<size_t>, set<vector<size_t>>> &equivalent_map,
                         map<vector<size_t>, int> &out_triplet_to_class)
{
  set<vector<size_t>> visited;
  vector<set<vector<size_t>>> classes;

  for (const auto &kv : equivalent_map)
  {
    const vector<size_t> &start = kv.first;
    if (visited.find(start) != visited.end())
      continue;

    set<vector<size_t>> cls;
    queue<vector<size_t>> q;
    q.push(start);

    while (!q.empty())
    {
      vector<size_t> cur = q.front();
      q.pop();
      if (visited.find(cur) != visited.end())
        continue;

      visited.insert(cur);
      cls.insert(cur);

      // forward neighbors
      auto it = equivalent_map.find(cur);
      if (it != equivalent_map.end())
      {
        for (const auto &nbr : it->second)
          if (visited.find(nbr) == visited.end())
            q.push(nbr);
      }

      // reverse neighbors: other keys that list cur as equivalent
      for (const auto &other_kv : equivalent_map)
      {
        const vector<size_t> &other_key = other_kv.first;
        if (visited.find(other_key) != visited.end())
          continue;
        const auto &other_set = other_kv.second;
        if (other_set.find(cur) != other_set.end())
          q.push(other_key);
      }
    }

    classes.push_back(cls);
  }

  // produce mapping triplet -> classId
  out_triplet_to_class.clear();
  for (size_t cid = 0; cid < classes.size(); ++cid)
  {
    for (const auto &t : classes[cid])
    {
      out_triplet_to_class[t] = static_cast<int>(cid);
    }
  }

  return classes;
}

// ---------- Main: compute equivalent triplets and classes ------------------
//
// Assumes: cfg and lce exist and the following functions/types exist:
//  - get_bcc_symmetry_operations(cfg, lce, uniqueR, uniqueT, symOps)
//  - FindClustersWithinAllowedSites(cfg, 3, 3, lce) returning objects with GetLatticeIdVector()
//  - cfg.GetRelativePositionOfLattice(id) -> Vector3d
//  - cfg.GetCentralAtomLatticeId()
//  - get_canonical_distance_order(triplet) -> vector<int>
//
void equivalent_clusters_triplets(Config &cfg, vector<size_t> lce, size_t latticeId)
{
  const double tolerance = 1e-6;

  // 1) Get symmetry operations
  unordered_set<Matrix3d, Matrix3dHash> uniqueR;
  unordered_set<Vector3d, Vector3dHash2> uniqueT;
  vector<pair<Matrix3d, Vector3d>> symOps;
  get_bcc_symmetry_operations(cfg, lce, uniqueR, uniqueT, symOps);

  // 2) Build id -> frac position map (and ensure centralId present in lce)
  unordered_map<size_t, Vector3d> idToPos;
  for (auto id : lce)
  {
    idToPos[id] = cfg.GetRelativePositionOfLattice(id);
  }

  // size_t centralId = cfg.GetCentralAtomLatticeId();
  if (idToPos.find(latticeId) == idToPos.end())
  {
    // add central if not present
    lce.push_back(latticeId);
    idToPos[latticeId] = cfg.GetRelativePositionOfLattice(latticeId);
  }

  // 3) Build all triplets that include centralId (canonicalized)
  auto allClusters = FindClustersWithinAllowedSites(cfg, 3, 3, lce);
  vector<vector<size_t>> triplets;
  for (const auto &c : allClusters)
  {
    auto ids = c.GetLatticeIdVector();
    if (ids.size() != 3)
      continue;
    if (find(ids.begin(), ids.end(), latticeId) == ids.end())
      continue;
    triplets.push_back(canonicalize_triplet(ids));
  }

  // Deduplicate triplets (since canonicalization can create duplicates)
  sort(triplets.begin(), triplets.end());
  triplets.erase(unique(triplets.begin(), triplets.end()), triplets.end());

  // 4) Build quantized pos->id map for fast matching
  unordered_map<string, size_t> posToId;
  for (const auto &p : idToPos)
  {
    posToId[frac_key(p.second)] = p.first;
  }

  auto get_canonical_distance_order = [&](const vector<size_t> &triplet) -> vector<int>
  {
    // triplet contains 3 IDs, centralId is fixed (should be in triplet)
    // Find the two non-central IDs
    vector<size_t> others;
    for (auto id : triplet)
      if (id != latticeId)
        others.push_back(id);

    // Get distance orders: central->others and others->others
    int d1 = cfg.GetDistanceOrder(latticeId, others[0]);
    int d2 = cfg.GetDistanceOrder(latticeId, others[1]);
    int d3 = cfg.GetDistanceOrder(others[0], others[1]);

    vector<int> distOrders = {d1, d2, d3};
    std::sort(distOrders.begin(), distOrders.end());
    return distOrders;
  };

  // 5) Optionally filter symmetry operations to those that fix the central position
  vector<pair<Matrix3d, Vector3d>> ops_fix_central;
  Vector3d centralPos = idToPos[latticeId];
  for (const auto &op : symOps)
  {
    Vector3d tp = op.first * centralPos + op.second;
    // wrap to [0,1)
    for (int k = 0; k < 3; ++k)
    {
      tp[k] = tp[k] - floor(tp[k]);
      if (tp[k] < 0)
        tp[k] += 1.0;
    }
    if ((tp - centralPos).norm() < tolerance)
      ops_fix_central.push_back(op);
  }
  if (ops_fix_central.empty())
    ops_fix_central = symOps; // fallback if none fixed

  // 6) For each triplet apply each op to the whole triplet and collect mapped triplets
  map<vector<size_t>, set<vector<size_t>>> equivalent_map;
  for (const auto &orig : triplets)
  {
    set<vector<size_t>> eqSet;
    eqSet.insert(orig); // include itself

    // compute original distance order (assumes get_canonical_distance_order works with canonical triplet)
    vector<int> origDist = get_canonical_distance_order(orig);

    for (const auto &op : ops_fix_central)
    {
      vector<size_t> mappedIds;
      mappedIds.reserve(3);
      bool all_found = true;

      for (auto id : orig)
      {
        Vector3d p = idToPos[id];
        Vector3d q = op.first * p + op.second;
        // wrap fractional to [0,1)
        for (int k = 0; k < 3; ++k)
        {
          q[k] = q[k] - floor(q[k]);
          if (q[k] < 0)
            q[k] += 1.0;
        }

        string key = frac_key(q);
        auto it = posToId.find(key);
        if (it != posToId.end())
        {
          mappedIds.push_back(it->second);
        }
        else
        {
          // fallback tolerant search if quantized key doesn't match
          bool found = false;
          for (const auto &pp : idToPos)
          {
            if ((pp.second - q).norm() < tolerance)
            {
              mappedIds.push_back(pp.first);
              found = true;
              break;
            }
          }
          if (!found)
          {
            all_found = false;
            break;
          }
        }
      } // end for each id in orig

      if (!all_found)
        continue;

      // canonicalize mapped triplet
      vector<size_t> mappedCanon = canonicalize_triplet(mappedIds);

      // ensure centralId is present (if that is your invariant)
      if (find(mappedCanon.begin(), mappedCanon.end(), latticeId) == mappedCanon.end())
        continue;

      // check distance-order equivalence (optional but preserves your earlier invariant)
      if (get_canonical_distance_order(mappedCanon) == origDist)
      {
        eqSet.insert(mappedCanon);
      }
    } // end for each op

    equivalent_map[orig] = eqSet;
  } // end for each triplet

  // 7) Form equivalence classes (connected components) and mapping triplet -> classId
  map<vector<size_t>, int> triplet_to_class;
  auto classes = form_equivalence_classes(equivalent_map, triplet_to_class);

  // 8) Print results (or return them depending on how you want to integrate)
  cout << "Found " << classes.size() << " equivalence classes.\n\n";
  for (size_t cid = 0; cid < classes.size(); ++cid)
  {
    cout << "Equivalence Class " << cid << ":\n";
    for (const auto &t : classes[cid])
    {
      cout << " Cluster: ";
      for (auto id : t)
        cout << id << " ";
      cout << "\n";
    }
    cout << "\n";
  }

  // Example: print triplet -> class id map
  cout << "Triplet -> ClassId map:\n";
  for (const auto &kv : triplet_to_class)
  {
    cout << "  [" << triplet_to_string(kv.first) << "] -> " << kv.second << "\n";
  }
}






vector<pair<Matrix3d, Vector3d>> GetSymmetryOperations(
    const Config &config,
    const unordered_set<size_t> &latticeIdVector,
    const double symprec,
    const map<size_t, size_t> &latticeIdToIndexMap)
{
  // BCC primitive cell
  double lattice[3][3] = {
      {1.0, 0.0, 0.0},
      {0.0, 1.0, 0.0},
      {0.0, 0.0, 1.0}};

  int numAtoms = latticeIdVector.size();

  double (*positions)[3] = new double[numAtoms][3];
  int *types = new int[numAtoms];

  int i = 0;
  for (const auto &latticeId : latticeIdVector)
  {
    VectorXd relativePosition = config.GetRelativePositionOfLattice(latticeId);

    positions[i][0] = relativePosition(0);
    positions[i][1] = relativePosition(1);
    positions[i][2] = relativePosition(2);

    types[i] = 1;
    i += 1;
  }

  SpglibDataset *dataset = spg_get_dataset(lattice, positions, types, numAtoms, symprec);

  if (!dataset)
  {
    cerr << "Failed to get spglib dataset." << endl;
    exit(1);
  }

  cout << "Space group: " << dataset->international_symbol
       << " (#" << dataset->spacegroup_number << ")" << endl;

  cout << "Point group: " << dataset->pointgroup_symbol << endl;

  vector<pair<Matrix3d, Vector3d>> symmetryOperations;
  for (int i = 0; i < dataset->n_operations; ++i)
  {
    Matrix3d rot;
    Vector3d trans;

    for (int r = 0; r < 3; ++r)
    {
      trans(r) = dataset->translations[i][r];
      for (int c = 0; c < 3; ++c)
      {
        rot(r, c) = dataset->rotations[i][r][c];
      }
    }
    pair<Matrix3d, Vector3d> symOpPair = {rot, trans};
    symmetryOperations.emplace_back(symOpPair);
  }

  // Print
  map<size_t, vector<size_t>> equivalentSitesMap;

  i = 0;
  for (const auto &latticeId : latticeIdVector)
  {
    size_t key = dataset->equivalent_atoms[i];

    if (equivalentSitesMap.find(key) != equivalentSitesMap.end())
    {
      equivalentSitesMap[key].push_back(latticeId); // use push_back instead of insert
    }
    else
    {
      equivalentSitesMap[key] = {latticeId};
    }
    // cout << value << " - " << key << endl;
    i += 1;
  }

  for (auto group : equivalentSitesMap)
  {
    cout << group.first << " : { ";
    for (auto id : group.second)
    {
      cout << id << " ";
    }
    cout << " } " << endl;

    if (latticeIdToIndexMap.size() != 0)
    {
      vector<size_t> indexVector;
      cout << "Index : ";
      for (auto id : group.second)
      {
        indexVector.emplace_back(latticeIdToIndexMap.at(id));
      }
      std::sort(indexVector.begin(), indexVector.end());
      print1DVector(indexVector);
    }
  }

  spg_free_dataset(dataset);
  delete[] positions;

  return symmetryOperations;
}

void GetEquivalentClusters(
    const Config &config,
    const unordered_set<size_t> &latticeIdVector,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const double symprec)
{

  auto symmetryOperations = GetSymmetryOperations(config,
                                                  latticeIdVector,
                                                  symprec);

  // Map latticeId to Position
  unordered_map<size_t, Vector3d> latticeIdToPositionMap;

  // Map Position to latticeId
  // String key inplace or positions
  unordered_map<Vector3d, size_t, Vector3dHash2, Vector3dEqual> positionToLatticeIdMap;

  for (auto latticeId : latticeIdVector)
  {
    Vector3d position = config.GetRelativePositionOfLattice(latticeId);
    latticeIdToPositionMap[latticeId] = position;
    positionToLatticeIdMap[position] = latticeId;
  }

  map<vector<size_t>, set<vector<size_t>>> equivalentMap;
  for (const auto &cluster : latticeClusterSet)
  {
    set<vector<size_t>> equivalentSet;
    set<LatticeCluster> equivalentLatticeClusters;

    auto originalCluster = cluster.GetLatticeIdVector();

    if (originalCluster.size() <= 1)
    {
      // Skip the singlets and empty cluster for now;
      continue;
    }

    // Iterate over the symmetry operations to get the equivalent clusters

    for (const auto &operation : symmetryOperations)
    {
      vector<size_t> transformedCluster;
      transformedCluster.reserve(originalCluster.size());
      bool allFound = true;

      // Tranform the cluster
      for (auto latticeId : originalCluster)
      {
        Vector3d position = latticeIdToPositionMap[latticeId];

        // operation = [rotation, translation]
        Vector3d transformedPosition = operation.first * position +
                                       operation.second;

        // wrap to [0, 1)
        for (int k = 0; k < 3; ++k)
        {
          transformedPosition[k] = transformedPosition[k] - floor(transformedPosition[k]);
          if (transformedPosition[k] < 0)
            transformedPosition[k] += 1.0;
        }

        // Get the lattice Id corresponding to the tranformed position
        auto it = positionToLatticeIdMap.find(transformedPosition);
        if (it != positionToLatticeIdMap.end())
        {
          transformedCluster.push_back(it->second);
        }
        // May be need to add other conditions
        else
        {
          allFound = false;
        }
      }

      if (!allFound)
        continue;

      auto latticeClusterType = IdentifyLatticeClusterType(config, transformedCluster);
      LatticeCluster tranformedLatticeCluster(latticeClusterType, transformedCluster);

      equivalentSet.insert(tranformedLatticeCluster.GetLatticeIdVector());
    }

    equivalentMap[originalCluster] = equivalentSet;
  }

  // Need to be rewritten

  map<vector<size_t>, int> triplet_to_class;
  auto classes = form_equivalence_classes(equivalentMap, triplet_to_class);

  // 8) Print results (or return them depending on how you want to integrate)
  cout << "Found " << classes.size() << " equivalence classes.\n\n";
  for (size_t cid = 0; cid < classes.size(); ++cid)
  {
    cout << "Equivalence Class " << cid << ":\n";
    for (const auto &t : classes[cid])
    {
      cout << " Cluster: ";
      for (auto id : t)
        cout << id << " ";
      auto clusterType = IdentifyLatticeClusterType(config, t);
      cout << clusterType << "\n";
    }

    cout << "\n";
  }
}

// Canonical LCE
// Choose one axis as representative
// Then for other axis apply symmetry operations to get the sorted atoms
// relative to the canonical representation

std::unordered_map<size_t, Eigen::RowVector3d> GetReference(const Config &config, size_t max_bond_order)
{

  // cout << "Reference :  0 - 25" << endl;
  // cout << config.GetRelativeDistanceVectorLattice(0, 25).transpose() << endl;

  pair<size_t, size_t> lattice_id_jump_pair = {0, 25};

  auto neighboring_lattice_ids = config.GetNeighboringLatticeIdSetOfPair(lattice_id_jump_pair,
                                                                         max_bond_order);

  size_t num_sites = neighboring_lattice_ids.size();
  Eigen::RowVector3d move_distance =
      Eigen::RowVector3d(0.5, 0.5, 0.5) - config.GetLatticePairCenter(lattice_id_jump_pair);

  std::unordered_map<size_t, Eigen::RowVector3d> lattice_id_hashmap;
  lattice_id_hashmap.reserve(num_sites);

  // Move lattice IDs to center
  for (const auto id : neighboring_lattice_ids)
  {
    Eigen::RowVector3d relative_position = config.GetRelativePositionOfLattice(id).transpose();
    relative_position += move_distance;
    relative_position -= relative_position.unaryExpr([](double x)
                                                     { return std::floor(x); });
    lattice_id_hashmap.emplace(id, relative_position);
  }

  return lattice_id_hashmap;
}

bool areVectorsClose(const Vector3d &a, const Vector3d &b, double tol = 1e-6)
{
  return (a - b).norm() < tol;
}

bool comparePositionMaps(const std::unordered_map<size_t, RowVector3d> &map1,
                         const std::unordered_map<size_t, RowVector3d> &map2,
                         double tol = 1e-6)
{
  if (map1.size() != map2.size())
    return false;

  // Extract positions into vectors
  std::vector<RowVector3d> vec1, vec2;
  for (const auto &[id, pos] : map1)
    vec1.push_back(pos);
  for (const auto &[id, pos] : map2)
    vec2.push_back(pos);

  // Custom comparator for sorting Vector3d lex order
  auto vectorLess = [](const RowVector3d &a, const RowVector3d &b)
  {
    if (std::abs(a.x() - b.x()) > 1e-8)
      return a.x() < b.x();
    if (std::abs(a.y() - b.y()) > 1e-8)
      return a.y() < b.y();
    return a.z() < b.z();
  };

  std::sort(vec1.begin(), vec1.end(), vectorLess);
  std::sort(vec2.begin(), vec2.end(), vectorLess);

  // Compare element-wise with tolerance
  for (size_t i = 0; i < vec1.size(); ++i)
  {
    if (!areVectorsClose(vec1[i], vec2[i], tol))
      return false;
  }

  return true;
}

void saveConfig(
    const Config &config,
    const unordered_map<size_t, Eigen::RowVector3d> &latticeIdPositionMap,
    const string &filename)
{

  size_t numAtoms = latticeIdPositionMap.size();

  Eigen::Matrix3Xd relativePositionMatrix(3, numAtoms);
  vector<Element> atomVector;
  atomVector.reserve(numAtoms);

  size_t colIndex = 0;
  for (const auto &entry : latticeIdPositionMap)
  {
    const Eigen::RowVector3d &pos = entry.second;
    // Eigen::RowVector3d is 1x3, but matrix expects 3x1 column, so transpose:
    relativePositionMatrix.col(colIndex) = pos.transpose();

    auto element = config.GetElementOfLattice(entry.first);

    atomVector.emplace_back(element);

    colIndex++;
  }

  auto cfg = Config(config.GetBasis(), relativePositionMatrix, atomVector);
  string filepath = "/home/pravendra3/Documents/LatticeMonteCarlo-eigen/bin/SymmetryOps/" + filename;
  Config::WriteConfig(filepath, cfg);
}

inline bool PositionCompareState(const std::pair<size_t, Eigen::RowVector3d> &lhs,
                                 const std::pair<size_t, Eigen::RowVector3d> &rhs)
{
  const auto &relative_position_lhs = lhs.second;
  const auto &relative_position_rhs = rhs.second;

  // Compare individual components (x, y, z)
  const double diff_x = relative_position_lhs[0] - relative_position_rhs[0];
  if (diff_x < -constants::kEpsilon)
  {
    return true;
  }
  if (diff_x > constants::kEpsilon)
  {
    return false;
  }

  const double diff_y = relative_position_lhs[1] - relative_position_rhs[1];
  if (diff_y < -constants::kEpsilon)
  {
    return true;
  }
  if (diff_y > constants::kEpsilon)
  {
    return false;
  }

  const double diff_z = relative_position_lhs[2] - relative_position_rhs[2];
  if (diff_z < -constants::kEpsilon)
  {
    return true;
  }
  if (diff_z > constants::kEpsilon)
  {
    return false;
  }

  return false;
}

vector<size_t> GetCanonicalSortedLatticeSites(
    const Config &config,
    const std::pair<size_t, size_t> &lattice_id_jump_pair,
    const size_t &max_bond_order)
{

  // One direction would be canonical direction 111
  // Using the reduced cords as reference just map the coordinates of other direction
  // along this 111 direction using symmetry operations

  // Canoncial Directon
  // 0 - 25 //  1 1 1 direction

  auto referenceLatticeIdHashmap = GetReference(
      config,
      max_bond_order);

  // for (auto entry : referenceLatticeIdHashmap)
  // {
  //   cout << entry.first << " : " << entry.second << endl;
  // }

  auto neighboring_lattice_ids = config.GetNeighboringLatticeIdSetOfPair(lattice_id_jump_pair,
                                                                         max_bond_order);

  size_t num_sites = neighboring_lattice_ids.size();
  Eigen::RowVector3d move_distance =
      Eigen::RowVector3d(0.5, 0.5, 0.5) - config.GetLatticePairCenter(lattice_id_jump_pair);

  std::unordered_map<size_t, Eigen::RowVector3d> lattice_id_hashmap;
  lattice_id_hashmap.reserve(num_sites);

  // Move lattice IDs to center
  for (const auto id : neighboring_lattice_ids)
  {
    Eigen::RowVector3d relative_position = config.GetRelativePositionOfLattice(id).transpose();
    relative_position += move_distance;
    relative_position -= relative_position.unaryExpr([](double x)
                                                     { return std::floor(x); });
    lattice_id_hashmap.emplace(id, relative_position);
  }

  // Symmetry operation for BCC
  auto latticeIdVectorSym = config.GetNeighborLatticeIdVectorOfLattice(0, 1);
  latticeIdVectorSym.emplace_back(0);

  unordered_set<size_t> latticeIdSet(latticeIdVectorSym.begin(), latticeIdVectorSym.end());
  auto symOps = GetSymmetryOperations(config, latticeIdSet);

  unordered_map<size_t, Eigen::RowVector3d> matchedTransformedPositionMap;

  for (auto &operation : symOps)
  {

    unordered_map<size_t, Eigen::RowVector3d> transformedPositionMap;

    for (auto entry : lattice_id_hashmap)
    {
      Vector3d position = entry.second;

      // operation = [rotation, translation]
      Vector3d transformedPosition = operation.first * position +
                                     operation.second;

      // wrap to [0, 1)
      for (int k = 0; k < 3; ++k)
      {
        transformedPosition[k] = transformedPosition[k] - floor(transformedPosition[k]);
        if (transformedPosition[k] < 0)
          transformedPosition[k] += 1.0;
      }

      transformedPositionMap[entry.first] = transformedPosition;
    }

    // Now compare the transformedPositionMap positions and referenceLatticeIdHashmap position

    auto isSame = comparePositionMaps(transformedPositionMap, referenceLatticeIdHashmap);
    if (isSame)
    {
      matchedTransformedPositionMap = transformedPositionMap;
      // cout << "Are Same" << endl;
      // cout << operation.first << endl;
      // cout << operation.second << endl;
      cout << endl;
      // for (auto entry: transformedPositionMap)
      // {
      //   cout << entry.first << " : " << entry.second << endl;
      // }
      break;
    }
  }

  // Save both the original and transformed one

  string filenameOriginal = "jumpPath_" + to_string(lattice_id_jump_pair.first) + "_" + to_string(lattice_id_jump_pair.second) + "_Original.cfg.gz";

  saveConfig(config, lattice_id_hashmap, filenameOriginal);
  string filenameTransformed = "jumpPath_" + to_string(lattice_id_jump_pair.first) + "_" + to_string(lattice_id_jump_pair.second) + "_Transformed.cfg.gz";

  saveConfig(config, matchedTransformedPositionMap, filenameTransformed);

  // sort the transformed positions

  // Convert unordered_map to vector for sorting
  std::vector<std::pair<size_t, Eigen::RowVector3d>>
      lattice_id_vector(matchedTransformedPositionMap.begin(), matchedTransformedPositionMap.end());

  // Sort the lattice vector based on PositionCompare
  std::sort(lattice_id_vector.begin(), lattice_id_vector.end(), PositionCompareState);

  // Extract and return only the lattice IDs
  std::vector<size_t> sorted_lattice_ids;
  sorted_lattice_ids.reserve(lattice_id_vector.size());
  for (const auto &pair : lattice_id_vector)
  {
    sorted_lattice_ids.push_back(pair.first);
    // std::cout << pair.first << " " << pair.second << std::endl;
  }

  return sorted_lattice_ids;
}
