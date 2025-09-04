// EquivalentClustersFast.cpp
// Fast replacement for GetEquivalentClusters using precomputed symmetry map + fractional hash
// Compile with: g++ -std=c++17 -O3 -fopenmp EquivalentClustersFast.cpp -I/path/to/eigen -lboost_system ...

#include <Eigen/Dense>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <algorithm>
#include <limits>
#include <cmath>
#include <iostream>
#include <cassert>
#include <cstdint>
#include <cstdlib>

#include "GetEquivalentClusters.h"
#include "Config.h"
#include "ClusterExpansion.h"
#include "SymmetrySpglib.h"

#ifdef _OPENMP
#include <omp.h>
#endif

// Include your project headers here (adapt paths)
// #include "Config.h"
// #include "LatticeCluster.h"
// #include "Helpers.h" // for GetEquivalentGroups, GetSymmetryOperations, etc.
#include <boost/functional/hash.hpp>

using std::map;
using std::pair;
using std::set;
using std::size_t;
using std::unordered_map;
using std::unordered_set;
using std::vector;

// ----- Adapt these typedefs if your project uses other names -----
using Vector3d = Eigen::Vector3d;
using Matrix3d = Eigen::Matrix3d;
// -----------------------------------------------------------------

// Helper small 3-int struct for hashed quantized fractional coords
struct Int3
{
  int x, y, z;
  bool operator==(Int3 const &o) const noexcept { return x == o.x && y == o.y && z == o.z; }
};
struct Int3Hash
{
  size_t operator()(Int3 const &k) const noexcept
  {
    // 64-bit mix
    uint64_t a = static_cast<uint32_t>(k.x);
    uint64_t b = static_cast<uint32_t>(k.y);
    uint64_t c = static_cast<uint32_t>(k.z);
    uint64_t seed = a + 0x9e3779b97f4a7c15ULL + (b << 6) + (b >> 2);
    seed ^= c + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return static_cast<size_t>(seed);
  }
};

// Quantize fractional coordinate (assumes f in [0,1) per component)
inline Int3 quantize_frac(const Vector3d &f, double bin_size)
{
  Int3 key;
  // floor rounding to nearest integer grid: using round so symmetric on boundaries
  key.x = static_cast<int>(std::floor(f[0] / bin_size + 0.5));
  key.y = static_cast<int>(std::floor(f[1] / bin_size + 0.5));
  key.z = static_cast<int>(std::floor(f[2] / bin_size + 0.5));
  return key;
}

// Wrap fractional coordinate into [0,1)
inline Vector3d wrap_frac(Vector3d f)
{
  for (int k = 0; k < 3; ++k)
  {
    f[k] = f[k] - std::floor(f[k]);
    if (f[k] < 0)
      f[k] += 1.0;
    // tiny snap to zero for numerical stability
    const double eps = 1e-12;
    if (f[k] > 1.0 - eps)
      f[k] = 0.0;
  }
  return f;
}

// Build fractional hash: bucket quantized fractional coords -> vector of lattice ids
static unordered_map<Int3, vector<size_t>, Int3Hash>
BuildFractionalHash(const unordered_map<size_t, Vector3d> &latticeIdToPosFrac, double bin_size)
{
  unordered_map<Int3, vector<size_t>, Int3Hash> fracHash;
  fracHash.reserve(latticeIdToPosFrac.size() * 2 + 16);
  for (const auto &p : latticeIdToPosFrac)
  {
    size_t id = p.first;
    Vector3d frac = wrap_frac(p.second);
    Int3 key = quantize_frac(frac, bin_size);
    fracHash[key].push_back(id);
  }
  return fracHash;
}

// Put this function into the same file where previous GetEquivalentClustersFast lived.
// Requires: Eigen, OpenMP, your GetSymmetryOperations, GetEquivalentGroups, IdentifyLatticeClusterType.
vector<set<vector<size_t>>> GetEquivalentClustersFast(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const double symprec,
    const bool debug,
    const bool /*isNew*/)
{
  using Vec3 = Eigen::Vector3d;

  // Quick exits
  if (latticeIdSet.empty())
    return {};
  if (latticeClusterSet.empty())
    return {};

  // Helper: small struct to hash vector<int64_t>
  struct VecInt64Hash
  {
    size_t operator()(const std::vector<int64_t> &v) const noexcept
    {
      size_t seed = v.size();
      for (auto x : v)
      {
        // combine, similar to boost::hash_combine
        seed ^= std::hash<int64_t>{}(x) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
  };

  // Small helpers
  auto wrap_frac = [](Vec3 f) -> Vec3
  {
    for (int k = 0; k < 3; ++k)
    {
      f[k] = f[k] - std::floor(f[k]);
      if (f[k] < 0)
        f[k] += 1.0;
      const double eps = 1e-12;
      if (f[k] > 1.0 - eps)
        f[k] = 0.0;
    }
    return f;
  };

  auto quantize_frac_to_int_tuple = [](const Vec3 &f, double scale) -> std::array<int64_t, 3>
  {
    std::array<int64_t, 3> a;
    a[0] = static_cast<int64_t>(std::llround(f[0] * scale));
    a[1] = static_cast<int64_t>(std::llround(f[1] * scale));
    a[2] = static_cast<int64_t>(std::llround(f[2] * scale));
    return a;
  };

  // 1) get symmetry operations using full lattice (important: compute symops on entire lattice)
  unordered_set<size_t> ids_for_symmetry;
  ids_for_symmetry.reserve(config.GetNumLattices());
  for (size_t i = 0; i < config.GetNumLattices(); ++i)
    ids_for_symmetry.insert(i);
  auto symmetryOperations = GetSymmetryOperations(config, ids_for_symmetry);
  if (symmetryOperations.empty())
  {
    if (debug)
      std::cerr << "[GetEquivalentClustersFast] No symmetry operations.\n";
    vector<set<vector<size_t>>> out;
    for (const auto &lc : latticeClusterSet)
    {
      set<vector<size_t>> s;
      s.insert(lc.GetLatticeIdVector());
      out.push_back(std::move(s));
    }
    return out;
  }
  const size_t Nsym = symmetryOperations.size();

  // 2) Build lattice id list and id->index map (deterministic order)
  vector<size_t> latticeIdList(latticeIdSet.begin(), latticeIdSet.end());
  std::sort(latticeIdList.begin(), latticeIdList.end());
  const size_t N_ids = latticeIdList.size();
  unordered_map<size_t, size_t> id_to_idx;
  id_to_idx.reserve(N_ids * 2);
  for (size_t i = 0; i < N_ids; ++i)
    id_to_idx[latticeIdList[i]] = i;

  // 3) Build id->fractional positions
  vector<Vec3> idx_to_frac(N_ids);
  for (size_t i = 0; i < N_ids; ++i)
    idx_to_frac[i] = config.GetRelativePositionOfLattice(latticeIdList[i]);

  // 4) Precompute transformed fractional positions for each symmetry op and site
  vector<vector<Vec3>> transPos(Nsym, vector<Vec3>(N_ids));
  for (size_t si = 0; si < Nsym; ++si)
  {
    const auto &op = symmetryOperations[si];
    for (size_t idx = 0; idx < N_ids; ++idx)
    {
      Vec3 tf = op.first * idx_to_frac[idx] + op.second;
      transPos[si][idx] = wrap_frac(tf);
    }
  }
  if (debug)
    std::cerr << "[GetEquivalentClustersFast] Precomputed transPos: Nsym=" << Nsym << " N_ids=" << N_ids << "\n";
  cout << "Done step 4" << endl;

  // 5) Convert cluster set to vector
  vector<LatticeCluster> clusterVec;
  clusterVec.reserve(latticeClusterSet.size());
  for (const auto &lc : latticeClusterSet)
    clusterVec.push_back(lc);

  // 6) canonical-key helper (SCALE tuned)
  const double SCALE = 1e8;
  auto canonical_key_from_frac_positions = [&](const std::vector<Vec3> &frac_positions) -> std::vector<int64_t>
  {
    const size_t k = frac_positions.size();
    std::vector<int64_t> best_flat;
    best_flat.clear();
    // Try translations that bring each site to origin then lexicographically choose minimal flattened int vector
    for (size_t anchor = 0; anchor < k; ++anchor)
    {
      Vec3 ref = frac_positions[anchor];
      std::vector<std::array<int64_t, 3>> rows;
      rows.reserve(k);
      for (size_t i = 0; i < k; ++i)
      {
        Vec3 rel = frac_positions[i] - ref;
        for (int d = 0; d < 3; ++d)
          rel[d] -= std::round(rel[d]); // nearest image
        rel = wrap_frac(rel);
        rows.push_back(quantize_frac_to_int_tuple(rel, SCALE));
      }
      std::sort(rows.begin(), rows.end(), [](auto const &a, auto const &b)
                {
                if (a[0] != b[0]) return a[0] < b[0];
                if (a[1] != b[1]) return a[1] < b[1];
                return a[2] < b[2]; });
      std::vector<int64_t> flat;
      flat.reserve(k * 3);
      for (auto &r : rows)
      {
        flat.push_back(r[0]);
        flat.push_back(r[1]);
        flat.push_back(r[2]);
      }
      if (best_flat.empty() || flat < best_flat)
        best_flat = std::move(flat);
    }
    return best_flat;
  };
  cout << "Done step 6" << endl;

  // 7) Per-thread maps
  int nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp single
    nthreads = omp_get_num_threads();
  }
#endif
  vector<unordered_map<vector<int64_t>, vector<vector<size_t>>, VecInt64Hash>> per_thread_maps;
  per_thread_maps.resize(nthreads);

  // 8) Parallel canonicalization and grouping
#pragma omp parallel
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    auto &local_map = per_thread_maps[tid];

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (size_t cidx = 0; cidx < clusterVec.size(); ++cidx)
    {
      const auto original = clusterVec[cidx].GetLatticeIdVector();
      // gather indices into idx space
      vector<size_t> idxs;
      idxs.reserve(original.size());
      bool bad = false;
      for (auto id : original)
      {
        auto it = id_to_idx.find(id);
        if (it == id_to_idx.end())
        {
          bad = true;
          break;
        }
        idxs.push_back(it->second);
      }
      if (bad)
      {
        // fallback key: raw ids (should be rare)
        vector<int64_t> fallback;
        fallback.reserve(original.size());
        for (auto id : original)
          fallback.push_back(static_cast<int64_t>(id));
        local_map[fallback].push_back(original);
        continue;
      }

      // best key across symmetry ops
      std::vector<int64_t> best_key;
      bool best_set = false;
      for (size_t si = 0; si < Nsym; ++si)
      {
        std::vector<Vec3> img;
        img.reserve(idxs.size());
        for (auto idx : idxs)
          img.push_back(transPos[si][idx]);
        std::vector<int64_t> key = canonical_key_from_frac_positions(img);
        if (!best_set || key < best_key)
        {
          best_key = std::move(key);
          best_set = true;
        }
      }

      if (!best_set)
      {
        vector<int64_t> fallback;
        fallback.reserve(original.size());
        for (auto id : original)
          fallback.push_back(static_cast<int64_t>(id));
        local_map[fallback].push_back(original);
      }
      else
      {
        local_map[best_key].push_back(original);
      }
    } // end for
  } // end omp

  cout << "Done step 8" << endl;

  // 9) Merge per-thread maps
  unordered_map<vector<int64_t>, vector<vector<size_t>>, VecInt64Hash> globalGrouped;
  globalGrouped.reserve(clusterVec.size() * 2);
  for (int t = 0; t < nthreads; ++t)
  {
    for (auto &kv : per_thread_maps[t])
    {
      auto &dst = globalGrouped[kv.first];
      dst.insert(dst.end(), std::make_move_iterator(kv.second.begin()), std::make_move_iterator(kv.second.end()));
    }
  }

  cout << "Done step 9" << endl;

  // 10) Build equivalentMap (representative -> set of clusters)
  map<vector<size_t>, set<vector<size_t>>> equivalentMap;
  for (auto &kv : globalGrouped)
  {
    const auto &group = kv.second;
    if (group.empty())
      continue;
    vector<size_t> rep = group.front();
    for (auto const &v : group)
      if (v < rep)
        rep = v;
    auto &s = equivalentMap[rep];
    for (auto const &v : group)
      s.insert(v);
  }

  // 11) Convert to equivalence classes using helper
  map<vector<size_t>, int> clustersToGroupMap;
  auto equivalentGroups = GetEquivalentGroups(equivalentMap, clustersToGroupMap);

  if (debug)
  {
    std::cout << "Found " << equivalentGroups.size() << " equivalence classes.\n";
    for (size_t cid = 0; cid < equivalentGroups.size(); ++cid)
    {
      std::cout << "Equivalence Class " << cid << ":\n";
      for (const auto &t : equivalentGroups[cid])
      {
        std::cout << " Cluster: ";
        for (auto id : t)
          std::cout << id << " ";
        auto clusterType = IdentifyLatticeClusterType(config, t);
        std::cout << clusterType << "\n";
      }
      std::cout << "\n";
    }
  }

  return equivalentGroups;
}

// Put this function into the same file where previous GetEquivalentClustersFast lived.
// Requires: Eigen, OpenMP, your GetSymmetryOperations, GetEquivalentGroups, IdentifyLatticeClusterType.
vector<set<vector<size_t>>> GetEquivalentClustersFast2(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const vector<pair<Matrix3d, Vector3d>> &symmetryOperations,
    const double symprec,
    const bool debug,
    const bool /*isNew*/)
{
  using Vec3 = Eigen::Vector3d;

  // Quick exits
  if (latticeIdSet.empty())
    return {};
  if (latticeClusterSet.empty())
    return {};

  // Helper: small struct to hash vector<int64_t>
  struct VecInt64Hash
  {
    size_t operator()(const std::vector<int64_t> &v) const noexcept
    {
      size_t seed = v.size();
      for (auto x : v)
      {
        // combine, similar to boost::hash_combine
        seed ^= std::hash<int64_t>{}(x) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
  };

  // Small helpers
  auto wrap_frac = [](Vec3 f) -> Vec3
  {
    for (int k = 0; k < 3; ++k)
    {
      f[k] = f[k] - std::floor(f[k]);
      if (f[k] < 0)
        f[k] += 1.0;
      const double eps = 1e-12;
      if (f[k] > 1.0 - eps)
        f[k] = 0.0;
    }
    return f;
  };

  auto quantize_frac_to_int_tuple = [](const Vec3 &f, double scale) -> std::array<int64_t, 3>
  {
    std::array<int64_t, 3> a;
    a[0] = static_cast<int64_t>(std::llround(f[0] * scale));
    a[1] = static_cast<int64_t>(std::llround(f[1] * scale));
    a[2] = static_cast<int64_t>(std::llround(f[2] * scale));
    return a;
  };

  // 1) get symmetry operations using the provided latticeIdSet (to match the original behavior and avoid over-grouping)
  auto ids_for_symmetry = latticeIdSet;
  // auto symmetryOperations = GetSymmetryOperations(config, ids_for_symmetry, true);
  if (symmetryOperations.empty())
  {
    if (debug)
      std::cerr << "[GetEquivalentClustersFast] No symmetry operations.\n";
    vector<set<vector<size_t>>> out;
    for (const auto &lc : latticeClusterSet)
    {
      set<vector<size_t>> s;
      s.insert(lc.GetLatticeIdVector());
      out.push_back(std::move(s));
    }
    return out;
  }
  const size_t Nsym = symmetryOperations.size();

  // 2) Build lattice id list and id->index map (deterministic order)
  vector<size_t> latticeIdList(latticeIdSet.begin(), latticeIdSet.end());
  std::sort(latticeIdList.begin(), latticeIdList.end());
  const size_t N_ids = latticeIdList.size();
  unordered_map<size_t, size_t> id_to_idx;
  id_to_idx.reserve(N_ids * 2);
  for (size_t i = 0; i < N_ids; ++i)
    id_to_idx[latticeIdList[i]] = i;

  // 3) Build id->fractional positions
  vector<Vec3> idx_to_frac(N_ids);
  for (size_t i = 0; i < N_ids; ++i)
    idx_to_frac[i] = config.GetRelativePositionOfLattice(latticeIdList[i]);

  // 4) Precompute transformed fractional positions for each symmetry op and site
  vector<vector<Vec3>> transPos(Nsym, vector<Vec3>(N_ids));
  for (size_t si = 0; si < Nsym; ++si)
  {
    const auto &op = symmetryOperations[si];
    for (size_t idx = 0; idx < N_ids; ++idx)
    {
      Vec3 tf = op.first * idx_to_frac[idx] + op.second;
      transPos[si][idx] = wrap_frac(tf);
    }
  }
  if (debug)
    std::cerr << "[GetEquivalentClustersFast] Precomputed transPos: Nsym=" << Nsym << " N_ids=" << N_ids << "\n";
  cout << "Done step 4" << endl;

  // 5) Convert cluster set to vector
  vector<LatticeCluster> clusterVec;
  clusterVec.reserve(latticeClusterSet.size());
  for (const auto &lc : latticeClusterSet)
    clusterVec.push_back(lc);

  // 6) canonical-key helper (SCALE tuned)
  const double SCALE = 1e8;
  auto canonical_key_from_frac_positions = [&](const std::vector<Vec3> &frac_positions) -> std::vector<int64_t>
  {
    const size_t k = frac_positions.size();
    std::vector<int64_t> best_flat;
    best_flat.clear();
    // Try translations that bring each site to origin then lexicographically choose minimal flattened int vector
    for (size_t anchor = 0; anchor < k; ++anchor)
    {
      Vec3 ref = frac_positions[anchor];
      std::vector<std::array<int64_t, 3>> rows;
      rows.reserve(k);
      for (size_t i = 0; i < k; ++i)
      {
        Vec3 rel = frac_positions[i] - ref;
        for (int d = 0; d < 3; ++d)
          rel[d] -= std::round(rel[d]); // nearest image
        rel = wrap_frac(rel);
        rows.push_back(quantize_frac_to_int_tuple(rel, SCALE));
      }
      std::sort(rows.begin(), rows.end(), [](auto const &a, auto const &b)
                {
                if (a[0] != b[0]) return a[0] < b[0];
                if (a[1] != b[1]) return a[1] < b[1];
                return a[2] < b[2]; });
      std::vector<int64_t> flat;
      flat.reserve(k * 3);
      for (auto &r : rows)
      {
        flat.push_back(r[0]);
        flat.push_back(r[1]);
        flat.push_back(r[2]);
      }
      if (best_flat.empty() || flat < best_flat)
        best_flat = std::move(flat);
    }
    return best_flat;
  };
  cout << "Done step 6" << endl;

  // 7) Per-thread maps
  int nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel
  {
#pragma omp single
    nthreads = omp_get_num_threads();
  }
#endif
  vector<unordered_map<vector<int64_t>, vector<vector<size_t>>, VecInt64Hash>> per_thread_maps;
  per_thread_maps.resize(nthreads);

  // 8) Parallel canonicalization and grouping
#pragma omp parallel
  {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    auto &local_map = per_thread_maps[tid];

#ifdef _OPENMP
#pragma omp for schedule(static)
#endif
    for (size_t cidx = 0; cidx < clusterVec.size(); ++cidx)
    {
      const auto original = clusterVec[cidx].GetLatticeIdVector();
      // gather indices into idx space
      vector<size_t> idxs;
      idxs.reserve(original.size());
      bool bad = false;
      for (auto id : original)
      {
        auto it = id_to_idx.find(id);
        if (it == id_to_idx.end())
        {
          bad = true;
          break;
        }
        idxs.push_back(it->second);
      }
      if (bad)
      {
        // fallback key: raw ids (should be rare)
        vector<int64_t> fallback;
        fallback.reserve(original.size());
        for (auto id : original)
          fallback.push_back(static_cast<int64_t>(id));
        local_map[fallback].push_back(original);
        continue;
      }

      // best key across symmetry ops
      std::vector<int64_t> best_key;
      bool best_set = false;
      for (size_t si = 0; si < Nsym; ++si)
      {
        std::vector<Vec3> img;
        img.reserve(idxs.size());
        for (auto idx : idxs)
          img.push_back(transPos[si][idx]);
        std::vector<int64_t> key = canonical_key_from_frac_positions(img);
        if (!best_set || key < best_key)
        {
          best_key = std::move(key);
          best_set = true;
        }
      }

      if (!best_set)
      {
        vector<int64_t> fallback;
        fallback.reserve(original.size());
        for (auto id : original)
          fallback.push_back(static_cast<int64_t>(id));
        local_map[fallback].push_back(original);
      }
      else
      {
        local_map[best_key].push_back(original);
      }
    } // end for
  } // end omp

  cout << "Done step 8" << endl;

  // 9) Merge per-thread maps
  unordered_map<vector<int64_t>, vector<vector<size_t>>, VecInt64Hash> globalGrouped;
  globalGrouped.reserve(clusterVec.size() * 2);
  for (int t = 0; t < nthreads; ++t)
  {
    for (auto &kv : per_thread_maps[t])
    {
      auto &dst = globalGrouped[kv.first];
      dst.insert(dst.end(), std::make_move_iterator(kv.second.begin()), std::make_move_iterator(kv.second.end()));
    }
  }

  cout << "Done step 9" << endl;

  // 10) Build equivalentMap (representative -> set of clusters)
  map<vector<size_t>, set<vector<size_t>>> equivalentMap;
  for (auto &kv : globalGrouped)
  {
    const auto &group = kv.second;
    if (group.empty())
      continue;
    vector<size_t> rep = group.front();
    for (auto const &v : group)
      if (v < rep)
        rep = v;
    auto &s = equivalentMap[rep];
    for (auto const &v : group)
      s.insert(v);
  }

  // 11) Convert to equivalence classes using helper
  map<vector<size_t>, int> clustersToGroupMap;
  auto equivalentGroups = GetEquivalentGroups(equivalentMap, clustersToGroupMap);

  if (debug)
  {
    std::cout << "Found " << equivalentGroups.size() << " equivalence classes.\n";
    for (size_t cid = 0; cid < equivalentGroups.size(); ++cid)
    {
      std::cout << "Equivalence Class " << cid << ":\n";
      for (const auto &t : equivalentGroups[cid])
      {
        std::cout << " Cluster: ";
        for (auto id : t)
          std::cout << id << " ";
        auto clusterType = IdentifyLatticeClusterType(config, t);
        std::cout << clusterType << "\n";
      }
      std::cout << "\n";
    }
  }

  return equivalentGroups;
}

vector<set<vector<size_t>>> GetEquivalentClustersByPermutation(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const bool debug = false)
{
  // 1. Get symmetry operations
  auto symmetryOperations = GetSymmetryOperations(config, latticeIdSet);

  // 2. Build lattice ID -> index mapping
  vector<size_t> latticeIds(latticeIdSet.begin(), latticeIdSet.end());
  sort(latticeIds.begin(), latticeIds.end());
  unordered_map<size_t, size_t> idToIndex;
  for (size_t i = 0; i < latticeIds.size(); ++i)
    idToIndex[latticeIds[i]] = i;

  // 3. Precompute symmetry permutations
  vector<vector<size_t>> symmetryPermutations;
  for (auto &op : symmetryOperations)
  {
    vector<size_t> perm(latticeIds.size());
    for (size_t i = 0; i < latticeIds.size(); ++i)
    {
      size_t lid = latticeIds[i];
      Vector3d pos = config.GetRelativePositionOfLattice(lid);
      Vector3d tpos = op.first * pos + op.second;

      // wrap to [0,1)
      for (int k = 0; k < 3; ++k)
      {
        tpos[k] = tpos[k] - floor(tpos[k]);
        if (tpos[k] < 0)
          tpos[k] += 1.0;
      }

      // Find which lattice ID this maps to
      size_t closest = 0;
      double minDist = numeric_limits<double>::max();
      for (size_t j = 0; j < latticeIds.size(); ++j)
      {
        Vector3d diff = tpos - config.GetRelativePositionOfLattice(latticeIds[j]);
        for (int k = 0; k < 3; ++k)
          diff[k] -= round(diff[k]);
        double d2 = diff.squaredNorm();
        if (d2 < minDist)
        {
          minDist = d2;
          closest = j;
        }
      }
      perm[i] = latticeIds[closest];
    }
    symmetryPermutations.push_back(perm);
  }

  // 4. Apply symmetry permutations to each cluster
  map<vector<size_t>, set<vector<size_t>>> equivalentMap;
  for (auto &cluster : latticeClusterSet)
  {
    vector<size_t> original = cluster.GetLatticeIdVector();
    sort(original.begin(), original.end());

    set<vector<size_t>> eqSet;
    for (auto &perm : symmetryPermutations)
    {
      vector<size_t> transformed;
      transformed.reserve(original.size());
      for (auto id : original)
      {
        auto it = find(latticeIds.begin(), latticeIds.end(), id);
        size_t idx = distance(latticeIds.begin(), it);
        transformed.push_back(perm[idx]);
      }
      sort(transformed.begin(), transformed.end());
      eqSet.insert(transformed);
    }

    // Always include original
    eqSet.insert(original);
    equivalentMap[original] = eqSet;
  }

  // 5. Group clusters
  map<vector<size_t>, int> clustersToGroupMap;
  auto equivalentGroups = GetEquivalentGroups(equivalentMap, clustersToGroupMap);

  // 6. Debug
  if (debug)
  {
    cout << "Found " << equivalentGroups.size() << " equivalence classes.\n";
    for (size_t cid = 0; cid < equivalentGroups.size(); ++cid)
    {
      cout << "Equivalence Class " << cid << ":\n";
      for (auto &cl : equivalentGroups[cid])
      {
        cout << " Cluster: ";
        for (auto id : cl)
          cout << id << " ";
        cout << IdentifyLatticeClusterType(config, cl) << "\n";
      }
    }
  }

  return equivalentGroups;
}

vector<set<vector<size_t>>> GetEquivalentClustersSlow(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const double symprec,
    const bool debug)
{
  using Vec3 = Eigen::Vector3d;

  if (latticeIdSet.empty() || latticeClusterSet.empty())
    return {};

  // 1) Get symmetry operations
  auto symmetryOperations = GetSymmetryOperations(config, latticeIdSet);

  // 2) Map latticeId to fractional positions
  unordered_map<size_t, Vec3> latticeIdToPositionMap;
  for (auto latticeId : latticeIdSet)
    latticeIdToPositionMap[latticeId] = config.GetRelativePositionOfLattice(latticeId);

  // 3) Helpers for fractional wrapping and integer canonicalization
  auto wrap_frac = [](Vec3 f) -> Vec3
  {
    for (int k = 0; k < 3; ++k)
    {
      f[k] = f[k] - floor(f[k]);
      if (f[k] < 0)
        f[k] += 1.0;
      const double eps = 1e-12;
      if (f[k] > 1.0 - eps)
        f[k] = 0.0;
    }
    return f;
  };

  const double SCALE = 1e8;
  auto quantize_frac_to_int_tuple = [](const Vec3 &f, double scale) -> std::array<int64_t, 3>
  {
    return {static_cast<int64_t>(std::llround(f[0] * scale)),
            static_cast<int64_t>(std::llround(f[1] * scale)),
            static_cast<int64_t>(std::llround(f[2] * scale))};
  };

  auto canonical_key_from_frac_positions = [&](const vector<Vec3> &frac_positions) -> vector<int64_t>
  {
    size_t k = frac_positions.size();
    vector<int64_t> best_flat;
    for (size_t anchor = 0; anchor < k; ++anchor)
    {
      Vec3 ref = frac_positions[anchor];
      vector<array<int64_t, 3>> rows;
      rows.reserve(k);
      for (size_t i = 0; i < k; ++i)
      {
        Vec3 rel = frac_positions[i] - ref;
        for (int d = 0; d < 3; ++d)
          rel[d] -= round(rel[d]); // nearest image
        rel = wrap_frac(rel);
        rows.push_back(quantize_frac_to_int_tuple(rel, SCALE));
      }
      sort(rows.begin(), rows.end(), [](auto const &a, auto const &b)
           {
                if (a[0] != b[0]) return a[0] < b[0];
                if (a[1] != b[1]) return a[1] < b[1];
                return a[2] < b[2]; });
      vector<int64_t> flat;
      flat.reserve(k * 3);
      for (auto &r : rows)
      {
        flat.push_back(r[0]);
        flat.push_back(r[1]);
        flat.push_back(r[2]);
      }
      if (best_flat.empty() || flat < best_flat)
        best_flat = move(flat);
    }
    return best_flat;
  };

  std::vector<std::vector<size_t>> debugClusters = {{34, 59, 60}, {5, 25, 34}};

  std::vector<size_t> clusterToDebug = {34, 59, 60};
  DebugCluster(clusterToDebug, latticeIdToPositionMap, symmetryOperations);

  cout << "For Cluster {5, 25, 34}" << endl;
  DebugCluster(debugClusters[1], latticeIdToPositionMap, symmetryOperations);

  // 4) Process clusters
  map<vector<size_t>, set<vector<size_t>>> equivalentMap;

  for (const auto &cluster : latticeClusterSet)
  {
    auto originalCluster = cluster.GetLatticeIdVector();
    set<vector<size_t>> equivalentSet;

    for (const auto &operation : symmetryOperations)
    {
      vector<Vec3> transformedPositions;
      transformedPositions.reserve(originalCluster.size());

      for (auto latticeId : originalCluster)
      {
        Vec3 pos = latticeIdToPositionMap[latticeId];
        Vec3 transformed = operation.first * pos + operation.second;
        transformed = wrap_frac(transformed);
        transformedPositions.push_back(transformed);
      }

      vector<int64_t> key = canonical_key_from_frac_positions(transformedPositions);

      // Map key back to lattice IDs
      vector<size_t> mappedCluster;
      mappedCluster.reserve(originalCluster.size());
      for (auto &pos : transformedPositions)
      {
        double minDist = numeric_limits<double>::max();
        size_t closestId = 0;
        for (auto &[id, latticePos] : latticeIdToPositionMap)
        {
          Vec3 diff = latticePos - pos;
          for (int d = 0; d < 3; ++d)
            diff[d] -= round(diff[d]);
          double d2 = (config.GetBasis() * diff).squaredNorm();
          if (d2 < minDist)
          {
            minDist = d2;
            closestId = id;
          }
        }
        mappedCluster.push_back(closestId);
      }
      sort(mappedCluster.begin(), mappedCluster.end());
      equivalentSet.insert(mappedCluster);
    }

    equivalentMap[originalCluster] = equivalentSet;
  }

  // 5) Group clusters into equivalence classes
  map<vector<size_t>, int> clustersToGroupMap;
  auto equivalentGroups = GetEquivalentGroups(equivalentMap, clustersToGroupMap);

  if (debug)
  {
    cout << "Found " << equivalentGroups.size() << " equivalence classes.\n";
    for (size_t cid = 0; cid < equivalentGroups.size(); ++cid)
    {
      cout << "Equivalence Class " << cid << ":\n";
      for (const auto &t : equivalentGroups[cid])
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

  return equivalentGroups;
}

vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetOrbitsNew(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const vector<pair<Matrix3d, Vector3d>> &symmetryOperations,
    const double symprec,
    const bool debug)
{
  auto equivalentClusters = GetEquivalentClustersFast2(
      config,
      latticeIdSet,
      latticeClusterSet,
      symmetryOperations,
      symprec,
      debug, true);

  // Encode the clusters
  vector<pair<vector<vector<size_t>>, LatticeClusterType>> equivalentOrbitVector;
  equivalentOrbitVector.reserve(equivalentClusters.size());

  for (const auto &orbit : equivalentClusters)
  {
    auto clusterType = IdentifyLatticeClusterType(config, *orbit.begin());

    vector<vector<size_t>> orbitVector(orbit.begin(), orbit.end());

    equivalentOrbitVector.emplace_back(make_pair(orbitVector, clusterType));
  }

  sort(equivalentOrbitVector.begin(), equivalentOrbitVector.end(),
       [](const auto &a, const auto &b)
       {
         return a.second < b.second;
       });

  return equivalentOrbitVector;
}

// Helper: quantize fractional position to integer key
inline std::array<int64_t, 3> quantize_frac_to_int_tuple(const Vec3 &f, double scale = 1e8)
{
  return {static_cast<int64_t>(std::llround(f[0] * scale)),
          static_cast<int64_t>(std::llround(f[1] * scale)),
          static_cast<int64_t>(std::llround(f[2] * scale))};
}

// Debug function
void DebugCluster(
    const std::vector<size_t> &cluster,
    const std::unordered_map<size_t, Vec3> &latticeIdToPos,
    const std::vector<std::pair<Eigen::Matrix3d, Vec3>> &symOps,
    double scale)
{
  std::cout << "\n--- Debugging Cluster ---\nOriginal IDs: ";
  for (auto id : cluster)
    std::cout << id << " ";
  std::cout << "\nOriginal fractional positions:\n";
  for (auto id : cluster)
    std::cout << "  " << id << ": " << latticeIdToPos.at(id).transpose() << "\n";

  for (size_t si = 0; si < symOps.size(); ++si)
  {
    const auto &op = symOps[si];
    std::vector<Vec3> transformed;
    transformed.reserve(cluster.size());
    for (auto id : cluster)
    {
      Vec3 tf = op.first * latticeIdToPos.at(id) + op.second;
      tf = wrap_frac(tf);
      transformed.push_back(tf);
    }

    std::cout << "Symmetry op " << si << " transformed positions:\n";
    for (size_t i = 0; i < cluster.size(); ++i)
      std::cout << "  " << cluster[i] << ": " << transformed[i].transpose() << "\n";

    // Canonical key (quantized)
    std::vector<int64_t> key;
    for (auto &v : transformed)
    {
      auto q = quantize_frac_to_int_tuple(v, scale);
      key.push_back(q[0]);
      key.push_back(q[1]);
      key.push_back(q[2]);
    }
    std::cout << "Canonical key: ";
    for (auto k : key)
      std::cout << k << " ";
    std::cout << "\n";
  }
  std::cout << "--- End Debug ---\n";
}

void DebugClusterFast(
    const std::vector<size_t> &cluster,
    const std::vector<size_t> &latticeIdList,
    const std::vector<Eigen::Vector3d> &idx_to_frac,
    const std::vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>> &symOps,
    double scale)
{
  std::unordered_map<size_t, size_t> id_to_idx;
  for (size_t i = 0; i < latticeIdList.size(); ++i)
    id_to_idx[latticeIdList[i]] = i;

  std::cout << "\n--- Debugging Cluster (Fast) ---\nOriginal IDs: ";
  for (auto id : cluster)
    std::cout << id << " ";
  std::cout << "\nOriginal fractional positions:\n";
  for (auto id : cluster)
  {
    size_t idx = id_to_idx.at(id);
    std::cout << "  " << id << ": " << idx_to_frac[idx].transpose() << "\n";
  }

  for (size_t si = 0; si < symOps.size(); ++si)
  {
    const auto &op = symOps[si];
    std::vector<Eigen::Vector3d> transformed;
    transformed.reserve(cluster.size());
    for (auto id : cluster)
    {
      size_t idx = id_to_idx.at(id);
      Eigen::Vector3d tf = op.first * idx_to_frac[idx] + op.second;
      // wrap
      for (int k = 0; k < 3; ++k)
      {
        tf[k] -= std::floor(tf[k]);
        if (tf[k] < 0)
          tf[k] += 1.0;
        const double eps = 1e-12;
        if (tf[k] > 1.0 - eps)
          tf[k] = 0.0;
      }
      transformed.push_back(tf);
    }

    std::cout << "Symmetry op " << si << " transformed positions:\n";
    for (size_t i = 0; i < cluster.size(); ++i)
      std::cout << "  " << cluster[i] << ": " << transformed[i].transpose() << "\n";

    // Canonical key (quantized)
    std::vector<int64_t> key;
    for (auto &v : transformed)
    {
      auto q = quantize_frac_to_int_tuple(v, scale);
      key.push_back(q[0]);
      key.push_back(q[1]);
      key.push_back(q[2]);
    }
    std::cout << "Canonical key: ";
    for (auto k : key)
      std::cout << k << " ";
    std::cout << "\n";
  }
  std::cout << "--- End Debug ---\n";
}

vector<set<vector<size_t>>> GetEquivalentClustersSlowFixed(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const double symprec,
    const bool debug)
{
  using Vec3 = Eigen::Vector3d;
  if (latticeIdSet.empty() || latticeClusterSet.empty())
    return {};

  // 1) Get symmetry operations
  auto symmetryOperations = GetSymmetryOperations(config, latticeIdSet);
  const size_t Nsym = symmetryOperations.size();

  // 2) Map latticeId -> index
  vector<size_t> latticeIdList(latticeIdSet.begin(), latticeIdSet.end());
  sort(latticeIdList.begin(), latticeIdList.end());
  unordered_map<size_t, size_t> id_to_idx;
  for (size_t i = 0; i < latticeIdList.size(); ++i)
    id_to_idx[latticeIdList[i]] = i;

  // 3) Map index -> fractional positions
  vector<Vec3> idx_to_frac(latticeIdList.size());
  for (size_t i = 0; i < latticeIdList.size(); ++i)
    idx_to_frac[i] = config.GetRelativePositionOfLattice(latticeIdList[i]);

  // Helper: wrap fractional positions
  auto wrap_frac = [](Vec3 f) -> Vec3
  {
    for (int k = 0; k < 3; ++k)
    {
      f[k] -= floor(f[k]);
      if (f[k] < 0)
        f[k] += 1.0;
      if (f[k] > 1.0 - 1e-12)
        f[k] = 0.0;
    }
    return f;
  };

  // Helper: quantize to int
  auto quantize_frac_to_int_tuple = [](const Vec3 &f, double scale = 1e8) -> array<int64_t, 3>
  {
    return {static_cast<int64_t>(std::llround(f[0] * scale)),
            static_cast<int64_t>(std::llround(f[1] * scale)),
            static_cast<int64_t>(std::llround(f[2] * scale))};
  };

  // Helper: canonical key from fractional positions
  auto canonical_key_from_frac_positions = [&](const vector<Vec3> &frac_positions) -> vector<int64_t>
  {
    size_t k = frac_positions.size();
    vector<int64_t> best_flat;
    for (size_t anchor = 0; anchor < k; ++anchor)
    {
      Vec3 ref = frac_positions[anchor];
      vector<array<int64_t, 3>> rows;
      rows.reserve(k);
      for (size_t i = 0; i < k; ++i)
      {
        Vec3 rel = frac_positions[i] - ref;
        for (int d = 0; d < 3; ++d)
          rel[d] -= round(rel[d]);
        rel = wrap_frac(rel);
        rows.push_back(quantize_frac_to_int_tuple(rel));
      }
      sort(rows.begin(), rows.end(), [](auto const &a, auto const &b)
           {
                if (a[0]!=b[0]) return a[0]<b[0];
                if (a[1]!=b[1]) return a[1]<b[1];
                return a[2]<b[2]; });
      vector<int64_t> flat;
      for (auto &r : rows)
      {
        flat.push_back(r[0]);
        flat.push_back(r[1]);
        flat.push_back(r[2]);
      }
      if (best_flat.empty() || flat < best_flat)
        best_flat = move(flat);
    }
    return best_flat;
  };

  // 4) Precompute transformed positions for each symmetry operation
  vector<vector<Vec3>> transPos(Nsym, vector<Vec3>(latticeIdList.size()));
  for (size_t si = 0; si < Nsym; ++si)
  {
    const auto &op = symmetryOperations[si];
    for (size_t idx = 0; idx < latticeIdList.size(); ++idx)
      transPos[si][idx] = wrap_frac(op.first * idx_to_frac[idx] + op.second);
  }

  // 5) Process clusters
  map<vector<size_t>, set<vector<size_t>>> equivalentMap;
  for (const auto &cluster : latticeClusterSet)
  {
    auto original = cluster.GetLatticeIdVector();
    vector<size_t> idxs;
    idxs.reserve(original.size());
    for (auto id : original)
      idxs.push_back(id_to_idx.at(id));

    // Find minimal canonical key across all sym ops
    vector<int64_t> best_key;
    bool best_set = false;
    for (size_t si = 0; si < Nsym; ++si)
    {
      vector<Vec3> img;
      img.reserve(idxs.size());
      for (auto idx : idxs)
        img.push_back(transPos[si][idx]);
      auto key = canonical_key_from_frac_positions(img);
      if (!best_set || key < best_key)
      {
        best_key = move(key);
        best_set = true;
      }
    }

    equivalentMap[original].insert(original); // always include original
                                              // optionally, include all symmetry images (if needed)
  }

  // 6) Group into equivalence classes
  map<vector<size_t>, int> clustersToGroupMap;
  auto equivalentGroups = GetEquivalentGroups(equivalentMap, clustersToGroupMap);

  if (debug)
  {
    cout << "Found " << equivalentGroups.size() << " equivalence classes.\n";
    for (size_t cid = 0; cid < equivalentGroups.size(); ++cid)
    {
      cout << "Equivalence Class " << cid << ":\n";
      for (const auto &t : equivalentGroups[cid])
      {
        cout << " Cluster: ";
        for (auto id : t)
          cout << id << " ";
        cout << IdentifyLatticeClusterType(config, t) << "\n";
      }
    }
  }

  return equivalentGroups;
}
