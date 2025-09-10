#include "GetEquivalentClusters.h"

static Vector3d ApplyPBCFractional(Vector3d f)
{
  for (int k = 0; k < 3; ++k)
  {
    f[k] -= floor(f[k]);
    if (f[k] < 0)
      f[k] += 1.0;
    if (f[k] > 1.0 - constants::kEpsilon)
      f[k] = 0.0;
  }
  return f;
}

static array<int64_t, 3> QuantizeFractionalToIntTuple(
    const Vector3d &f)
{
  array<int64_t, 3> a;
  for (int i = 0; i < 3; ++i)
  {
    a[i] = static_cast<int64_t>(llround(f[i] * constants::K_FRACTIONAL_TO_INT_SCALE));
  }
  return a;
}

static vector<int64_t> CanonicalKeyFromFractionalPosition(
    const vector<Vector3d> &frac_positions)
{
  const size_t k = frac_positions.size();
  vector<int64_t> best_flat;
  best_flat.clear();

  // Try translations that bring each site to origin then lexicographically choose minimal flattened int vector
  for (size_t anchor = 0; anchor < k; ++anchor)
  {
    Vector3d ref = frac_positions[anchor];
    vector<array<int64_t, 3>> rows;
    rows.reserve(k);

    for (size_t i = 0; i < k; ++i)
    {
      Vector3d rel = frac_positions[i] - ref;
      for (int d = 0; d < 3; ++d)
        rel[d] -= round(rel[d]); // nearest image

      rel = ApplyPBCFractional(rel);
      rows.push_back(QuantizeFractionalToIntTuple(rel));
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

// Helper function to get the equivalent clusters
static vector<set<vector<size_t>>> GetEquivalentGroups(
    const map<vector<size_t>, set<vector<size_t>>> &equivalentMap,
    map<vector<size_t>, int> &clustersToGroupMap)
{
  set<vector<size_t>> visited;
  vector<set<vector<size_t>>> groups;

  for (const auto &kv : equivalentMap)
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
      auto it = equivalentMap.find(cur);
      if (it != equivalentMap.end())
      {
        for (const auto &nbr : it->second)
          if (visited.find(nbr) == visited.end())
            q.push(nbr);
      }

      // reverse neighbors: other keys that list cur as equivalent
      for (const auto &other_kv : equivalentMap)
      {
        const vector<size_t> &other_key = other_kv.first;
        if (visited.find(other_key) != visited.end())
          continue;
        const auto &other_set = other_kv.second;
        if (other_set.find(cur) != other_set.end())
          q.push(other_key);
      }
    }

    groups.push_back(cls);
  }

  // produce mapping triplet -> classId
  clustersToGroupMap.clear();
  for (size_t cid = 0; cid < groups.size(); ++cid)
  {
    for (const auto &t : groups[cid])
    {
      clustersToGroupMap[t] = static_cast<int>(cid);
    }
  }

  return groups;
}

// Removes the specified lattice IDs from all clusters in a set of orbits
static vector<set<vector<size_t>>> GetLocalClustersExcludingSites(
    const vector<set<vector<size_t>>> &equivalentClusters,
    const vector<size_t> &excludedLatticeIds)
{
  vector<set<vector<size_t>>> localClusters;
  localClusters.reserve(equivalentClusters.size());

  for (const auto &orbit : equivalentClusters)
  {
    set<vector<size_t>> newOrbit;
    for (auto cluster : orbit)
    {
      for (auto id : excludedLatticeIds)
      {
        cluster.erase(remove(cluster.begin(), cluster.end(), id), cluster.end());
      }
      newOrbit.insert(cluster);
    }
    localClusters.emplace_back(newOrbit);
  }
  return localClusters;
}

vector<set<vector<size_t>>> GetEquivalentClusters(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const bool debug,
    const double symprec)
{
  // Quick exits and validation
  if (latticeIdSet.empty())
    throw runtime_error("Error in `GetEquivalentClusters`: latticeIdSet is empty.");

  if (latticeClusterSet.empty())
    throw runtime_error("Error in `GetEquivalentClusters`: latticeClusterSet is empty.");

  auto symmetryOperations = GetSymmetryOperations(
      config,
      latticeIdSet,
      debug);

  if (symmetryOperations.empty())
    throw runtime_error("Error in `GetEquivalentClusters`: symmetryOperations is empty.");

  // Number of symmetry operations
  const size_t numSymmetryOps = symmetryOperations.size();

  // Build lattice id list and id->index map (deterministic order)
  vector<size_t> latticeIdList(latticeIdSet.begin(), latticeIdSet.end());
  sort(latticeIdList.begin(), latticeIdList.end());

  const size_t N_ids = latticeIdList.size();

  // Lattice Id to Index Map
  unordered_map<size_t, size_t> id_to_idx;
  id_to_idx.reserve(N_ids * 2);
  for (size_t i = 0; i < N_ids; ++i)
    id_to_idx[latticeIdList[i]] = i;

  // Index to Fractional Position Map
  vector<Vector3d> idx_to_frac(N_ids);
  for (size_t i = 0; i < N_ids; ++i)
    idx_to_frac[i] = config.GetRelativePositionOfLattice(latticeIdList[i]);

  // Precompute transformed fractional positions for each symmetry op and site
  vector<vector<Vector3d>> transformedPositions(numSymmetryOps, vector<Vector3d>(N_ids));

  for (size_t si = 0; si < numSymmetryOps; ++si)
  {
    const auto &op = symmetryOperations[si];
    for (size_t idx = 0; idx < N_ids; ++idx)
    {
      Vector3d tf = op.first * idx_to_frac[idx] + op.second;
      transformedPositions[si][idx] = ApplyPBCFractional(tf);
    }
  }

  // Convert cluster set to vector
  vector<LatticeCluster> latticeClustersVector;

  latticeClustersVector.reserve(latticeClusterSet.size());
  for (const auto &latticeCluster : latticeClusterSet)
    latticeClustersVector.push_back(latticeCluster);

  // Canonicalization & Grouping
  unordered_map<vector<int64_t>, vector<vector<size_t>>, VecInt64Hash> groupedClusters;
  groupedClusters.reserve(latticeClustersVector.size() * 2);

  for (size_t cidx = 0; cidx < latticeClustersVector.size(); ++cidx)
  {
    const auto &original = latticeClustersVector[cidx].GetLatticeIdVector();

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
      vector<int64_t> fallback;
      fallback.reserve(original.size());
      for (auto id : original)
        fallback.push_back(static_cast<int64_t>(id));
      groupedClusters[fallback].push_back(original);
      continue;
    }

    // best key across symmetry ops
    vector<int64_t> best_key;
    bool best_set = false;
    for (size_t si = 0; si < numSymmetryOps; ++si)
    {
      vector<Vector3d> img;
      img.reserve(idxs.size());
      for (auto idx : idxs)
        img.push_back(transformedPositions[si][idx]);
      vector<int64_t> key = CanonicalKeyFromFractionalPosition(img);
      if (!best_set || key < best_key)
      {
        best_key = move(key);
        best_set = true;
      }
    }

    if (!best_set)
    {
      vector<int64_t> fallback;
      fallback.reserve(original.size());
      for (auto id : original)
        fallback.push_back(static_cast<int64_t>(id));
      groupedClusters[fallback].push_back(original);
    }
    else
    {
      groupedClusters[best_key].push_back(original);
    }
  }

  // Equivalent Map
  map<vector<size_t>, set<vector<size_t>>> equivalentMap;
  for (auto &kv : groupedClusters)
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

  // Convert to equivalence classes using GetEquivalentGroups
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

vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetEncodedOrbits(
    const Config &config,
    const vector<size_t> &sortedNNLatticeIdVector,
    const vector<set<vector<size_t>>> &equivalentClusters,
    const bool debug)
{
  map<size_t, size_t> latticeIdToIndexMap;

  for (size_t i = 0; i < sortedNNLatticeIdVector.size(); i++)
  {
    latticeIdToIndexMap[sortedNNLatticeIdVector[i]] = i;
  }

  // Encode the clusters
  vector<pair<vector<vector<size_t>>, LatticeClusterType>> encodedEquivalentClusters;
  encodedEquivalentClusters.reserve(equivalentClusters.size());

  for (const auto &orbits : equivalentClusters)
  {
    vector<vector<size_t>> encodedOrbits;
    LatticeClusterType clusterType;
    for (const auto &latticeCluster : orbits)
    {
      clusterType = IdentifyLatticeClusterType(config, latticeCluster);
      if (latticeCluster.size() == 0)
      {
        encodedOrbits.emplace_back(latticeCluster);
      }
      else
      {
        vector<size_t> encodedCluster;
        encodedCluster.reserve(latticeCluster.size());
        for (const auto latticeId : latticeCluster)
        {
          encodedCluster.emplace_back(latticeIdToIndexMap.at(latticeId));
        }
        encodedOrbits.emplace_back(encodedCluster);
      }
    }
    encodedEquivalentClusters.emplace_back(make_pair(encodedOrbits, clusterType));
  }

  sort(encodedEquivalentClusters.begin(), encodedEquivalentClusters.end(),
       [](const auto &a, const auto &b)
       {
         return a.second < b.second;
       });

  bool emptyOrbitKept = false;

  // This part takes care that only one orbit with empty cluster will be there 
  // so as to avoid redundancy in intercept term
  encodedEquivalentClusters.erase(
      std::remove_if(
          encodedEquivalentClusters.begin(),
          encodedEquivalentClusters.end(),
          [&](const auto &orbit)
          {
            // Check if orbit is entirely empty clusters
            bool allEmpty = !orbit.first.empty() &&
                            std::all_of(orbit.first.begin(), orbit.first.end(),
                                        [](const auto &cluster)
                                        {
                                          return cluster.empty();
                                        });

            if (allEmpty)
            {
              if (!emptyOrbitKept)
              {
                emptyOrbitKept = true; // keep the first one
                return false;          // don’t remove
              }
              return true; // remove subsequent ones
            }
            return false; // keep normal orbits
          }),
      encodedEquivalentClusters.end());

  if (debug)
  {
    cout << "Encoded Equivalent Clusters : " << endl;
    int orbitIndex = 0;
    for (const auto &pair : encodedEquivalentClusters)
    {
      const auto &encodedOrbits = pair.first;
      const auto &clusterType = pair.second;

      cout << "Cluster Type: " << clusterType << "\n";
      cout << "Orbits: " << orbitIndex << endl;

      for (const auto &orbit : encodedOrbits)
      {
        cout << "  [ ";
        for (const auto &id : orbit)
          cout << id << " ";
        cout << "]\n";
      }

      cout << "------------------------\n";
      orbitIndex++;
    }

    unordered_set<size_t> nnLatticeIdSet(
        sortedNNLatticeIdVector.begin(),
        sortedNNLatticeIdVector.end());
    GetSymmetryOperations(config, nnLatticeIdSet, debug, latticeIdToIndexMap);
  }

  return encodedEquivalentClusters;
}

vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetLocalEncodedOrbitsForSite(
    const Config &config,
    const size_t &maxBondOrder,
    const size_t &maxClusterSize,
    const bool debug,
    const double symprec)
{
  size_t latticeId = 0;

  // Canonical Sorting
  const auto sortedNNLatticeIdVector = GetCanonicalSortedSitesForSite(config, latticeId, maxBondOrder);

  unordered_set<size_t> nnLatticeIdSet(sortedNNLatticeIdVector.begin(), sortedNNLatticeIdVector.end());
  nnLatticeIdSet.insert(latticeId);

  // Symmetry Operations
  auto symmetryOperation = GetSymmetryOperations(config, nnLatticeIdSet, debug, {}, symprec);

  size_t newClusterSize = maxClusterSize + 1;

  // Lattice clusters centered around a site
  auto latticeClustersSet = FindAllLatticeClusters(
      config,
      newClusterSize,
      maxBondOrder,
      {latticeId});

  // Equivalent Clusters including the lattice site
  auto equivalentClusters = GetEquivalentClusters(
      config,
      nnLatticeIdSet,
      latticeClustersSet, symprec, debug);

  // Remove the central site from the clusters to get the local clusters
  // which defines the local environment around a latticeId
  auto localEquivalentClusters = GetLocalClustersExcludingSites(
      equivalentClusters,
      {latticeId});

  // Encode the local equivalent cluster using sortedNNLatticeIdVector
  // which only contains the neighbours of the latticeId

  return GetEncodedOrbits(
      config,
      sortedNNLatticeIdVector,
      localEquivalentClusters,
      debug);
}

vector<pair<vector<vector<size_t>>, LatticeClusterType>> GetLocalEncodedOrbitsForPair(
    const Config &config,
    const size_t &maxBondOrder,
    const size_t &maxClusterSize,
    const Vector3d &referenceJumpDirection,
    const bool debug,
    const double symprec)
{
  pair<size_t, size_t> jumpPair = {0, config.GetNeighborLatticeIdVectorOfLattice(0, 1)[0]};

  auto nnLatticeSites = config.GetNeighboringLatticeIdSetOfPair(jumpPair, maxBondOrder);
  nnLatticeSites.insert(jumpPair.first);
  nnLatticeSites.insert(jumpPair.second);

  auto symmetryOperation = GetSymmetryOperations(config, nnLatticeSites);

  // Will be used to get the canonical sorted lattice sites
  auto canonicalReferenceMap = GetCenteredNeighborsAlongJumpDirection(
      config,
      maxBondOrder,
      referenceJumpDirection);

  // canonically sorted lattice id based on the canonicalReferenceMap
  const auto canonicalSortedLatticeSites = GetCanonicalSortedSitesForPair(
      config,
      jumpPair,
      maxBondOrder,
      canonicalReferenceMap,
      symmetryOperation);

  size_t newClusterSize = maxClusterSize + 1;

  // Lattice clusters centered around a sites in the pair
  auto latticeClustersSet = FindAllLatticeClusters(
      config,
      newClusterSize,
      maxBondOrder,
      {jumpPair.first, jumpPair.second});

  // Equivalent cluster including the lattice sites in the pair
  auto equivalentClusters = GetEquivalentClusters(
      config,
      nnLatticeSites,
      latticeClustersSet, symprec, debug);

  // Remove the lattice sites in the pair from the clusters to get the
  // local clusters which defines the local environment around the jump pair
  auto localEquivalentClusters = GetLocalClustersExcludingSites(
      equivalentClusters,
      {jumpPair.first, jumpPair.second});

  // Encode the local equivalent cluster using canonicalSortedLatticeSites
  // which only contains the neighbours of the latticeId jump pair

  return GetEncodedOrbits(
      config,
      canonicalSortedLatticeSites,
      localEquivalentClusters,
      debug);
}

/// Original : Need to be removed
/*

vector<set<vector<size_t>>> GetEquivalentClusters(
    const Config &config,
    const unordered_set<size_t> &latticeIdSet,
    const unordered_set<LatticeCluster, boost::hash<LatticeCluster>> &latticeClusterSet,
    const bool debug,
    const double symprec)
{

  auto symmetryOperations = GetSymmetryOperations(
      config,
      latticeIdSet,
      debug);

  // Quick exits
  if (latticeIdSet.empty())
    return {};
  if (latticeClusterSet.empty())
    return {};

  // Helper: small struct to hash vector<int64_t>
  struct VecInt64Hash
  {
    size_t operator()(const vector<int64_t> &v) const noexcept
    {
      size_t seed = v.size();
      for (auto x : v)
      {
        // combine, similar to boost::hash_combine
        seed ^= hash<int64_t>{}(x) + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
      }
      return seed;
    }
  };

  // Small helpers
  auto ApplyPBCFractional = [](Vector3d f) -> Vector3d
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

  auto QuantizeFractionalToIntTuple = [](const Vector3d &f, double scale) -> array<int64_t, 3>
  {
    array<int64_t, 3> a;
    a[0] = static_cast<int64_t>(llround(f[0] * scale));
    a[1] = static_cast<int64_t>(llround(f[1] * scale));
    a[2] = static_cast<int64_t>(llround(f[2] * scale));
    return a;
  };

  // 1) get symmetry operations using the provided latticeIdSet (to match the original behavior and avoid over-grouping)
  auto ids_for_symmetry = latticeIdSet;
  // auto symmetryOperations = GetSymmetryOperations(config, ids_for_symmetry, true);
  if (symmetryOperations.empty())
  {
    if (debug)
      cerr << "[GetEquivalentClustersFast] No symmetry operations.\n";
    vector<set<vector<size_t>>> out;
    for (const auto &lc : latticeClusterSet)
    {
      set<vector<size_t>> s;
      s.insert(lc.GetLatticeIdVector());
      out.push_back(move(s));
    }
    return out;
  }
  const size_t numSymmetryOps = symmetryOperations.size();

  // 2) Build lattice id list and id->index map (deterministic order)
  vector<size_t> latticeIdList(latticeIdSet.begin(), latticeIdSet.end());
  sort(latticeIdList.begin(), latticeIdList.end());
  const size_t N_ids = latticeIdList.size();
  unordered_map<size_t, size_t> id_to_idx;
  id_to_idx.reserve(N_ids * 2);
  for (size_t i = 0; i < N_ids; ++i)
    id_to_idx[latticeIdList[i]] = i;

  // 3) Build id->fractional positions
  vector<Vector3d> idx_to_frac(N_ids);
  for (size_t i = 0; i < N_ids; ++i)
    idx_to_frac[i] = config.GetRelativePositionOfLattice(latticeIdList[i]);

  // 4) Precompute transformed fractional positions for each symmetry op and site
  vector<vector<Vector3d>> transformedPositions(numSymmetryOps, vector<Vector3d>(N_ids));
  for (size_t si = 0; si < numSymmetryOps; ++si)
  {
    const auto &op = symmetryOperations[si];
    for (size_t idx = 0; idx < N_ids; ++idx)
    {
      Vector3d tf = op.first * idx_to_frac[idx] + op.second;
      transformedPositions[si][idx] = ApplyPBCFractional(tf);
    }
  }
  if (debug)
    cerr << "[GetEquivalentClustersFast] Precomputed transformedPositions: numSymmetryOps=" << numSymmetryOps << " N_ids=" << N_ids << "\n";
  cout << "Done step 4" << endl;

  // 5) Convert cluster set to vector
  vector<LatticeCluster> latticeClustersVector;
  latticeClustersVector.reserve(latticeClusterSet.size());
  for (const auto &lc : latticeClusterSet)
    latticeClustersVector.push_back(lc);

  // 6) canonical-key helper (SCALE tuned)
  const double SCALE = 1e8;
  auto CanonicalKeyFromFractionalPosition = [&](const vector<Vector3d> &frac_positions) -> vector<int64_t>
  {
    const size_t k = frac_positions.size();
    vector<int64_t> best_flat;
    best_flat.clear();
    // Try translations that bring each site to origin then lexicographically choose minimal flattened int vector
    for (size_t anchor = 0; anchor < k; ++anchor)
    {
      Vector3d ref = frac_positions[anchor];
      vector<array<int64_t, 3>> rows;
      rows.reserve(k);
      for (size_t i = 0; i < k; ++i)
      {
        Vector3d rel = frac_positions[i] - ref;
        for (int d = 0; d < 3; ++d)
          rel[d] -= round(rel[d]); // nearest image
        rel = ApplyPBCFractional(rel);
        rows.push_back(QuantizeFractionalToIntTuple(rel, SCALE));
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
    for (size_t cidx = 0; cidx < latticeClustersVector.size(); ++cidx)
    {
      const auto original = latticeClustersVector[cidx].GetLatticeIdVector();
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
      vector<int64_t> best_key;
      bool best_set = false;
      for (size_t si = 0; si < numSymmetryOps; ++si)
      {
        vector<Vector3d> img;
        img.reserve(idxs.size());
        for (auto idx : idxs)
          img.push_back(transformedPositions[si][idx]);
        vector<int64_t> key = CanonicalKeyFromFractionalPosition(img);
        if (!best_set || key < best_key)
        {
          best_key = move(key);
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
  globalGrouped.reserve(latticeClustersVector.size() * 2);
  for (int t = 0; t < nthreads; ++t)
  {
    for (auto &kv : per_thread_maps[t])
    {
      auto &dst = globalGrouped[kv.first];
      dst.insert(dst.end(), make_move_iterator(kv.second.begin()), make_move_iterator(kv.second.end()));
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

*/