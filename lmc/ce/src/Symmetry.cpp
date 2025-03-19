#include "Symmetry.h"

inline bool PositionCompareState(
    const pair<size_t, RowVector3d> &lhs,
    const pair<size_t, RowVector3d> &rhs)
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

void RotateLatticeVector(
    unordered_map<size_t, RowVector3d> &lattice_id_hashmap,
    const Matrix3d &rotation_matrix)
{
  const RowVector3d
      move_distance_after_rotation = RowVector3d(0.5, 0.5, 0.5) -
                                     (RowVector3d(0.5, 0.5, 0.5) * rotation_matrix);

  for (auto &lattice : lattice_id_hashmap)
  {

    auto relative_position = lattice.second;
    // rotate
    relative_position = relative_position * rotation_matrix;
    // move to new center
    relative_position += move_distance_after_rotation;
    relative_position -= relative_position.unaryExpr([](double x)
                                                     { return floor(x); });

    lattice.second = relative_position;
  }
}

vector<size_t> GetSortedLatticeVectorStateOfPair(
    const Config &config,
    const pair<size_t, size_t> &lattice_id_jump_pair,
    const size_t &max_bond_order)
{
  auto neighboring_lattice_ids =
      config.GetNeighboringLatticeIdSetOfPair(lattice_id_jump_pair,
                                              max_bond_order);

  size_t num_sites = neighboring_lattice_ids.size();

  RowVector3d move_distance = RowVector3d(0.5, 0.5, 0.5) -
                              GetLatticePairCenter(config, lattice_id_jump_pair);

  unordered_map<size_t, RowVector3d> lattice_id_hashmap;
  lattice_id_hashmap.reserve(num_sites);

  // Move lattice IDs to center
  for (const auto id : neighboring_lattice_ids)
  {
    RowVector3d relative_position =
        config.GetRelativePositionOfLattice(id).transpose();

    relative_position += move_distance;
    relative_position -= relative_position.unaryExpr([](double x)
                                                     { return floor(x); });
    lattice_id_hashmap.emplace(id, relative_position);
  }

  // Rotate lattice vectors
  auto rotationMatrix = GetLatticePairRotationMatrix(config,
                                                     lattice_id_jump_pair);
  RotateLatticeVector(lattice_id_hashmap, rotationMatrix);

  // Convert unordered_map to vector for sorting
  vector<pair<size_t, RowVector3d>>
      lattice_id_vector(lattice_id_hashmap.begin(), lattice_id_hashmap.end());

  // Sort the lattice vector based on PositionCompareMMM
  sort(lattice_id_vector.begin(), lattice_id_vector.end(), PositionCompareState);

  // Extract and return only the lattice IDs
  vector<size_t> sorted_lattice_ids;
  sorted_lattice_ids.reserve(lattice_id_vector.size());
  for (const auto &pair : lattice_id_vector)
  {
    sorted_lattice_ids.push_back(pair.first);
  }

  return sorted_lattice_ids;
}

RowVector3d GetLatticePairCenter(
    const Config &config,
    const pair<size_t, size_t> &lattice_id_jump_pair)
{
  // Get relative positions of the lattice IDs
  Vector3d first_relative =
      config.GetRelativePositionOfLattice(lattice_id_jump_pair.first);

  Vector3d second_relative =
      config.GetRelativePositionOfLattice(lattice_id_jump_pair.second);

  Vector3d center_position;
  for (int kDim = 0; kDim < 3; ++kDim)
  {
    double distance = first_relative[kDim] - second_relative[kDim];
    int period = static_cast<int>(distance / 0.5);

    // Adjust the positions once based on the period
    if (period != 0)
    {
      first_relative[kDim] -= period;
    }

    center_position[kDim] = 0.5 * (first_relative[kDim] + second_relative[kDim]);
  }

  return center_position.transpose();
}

Matrix3d GetLatticePairRotationMatrix(
    const Config &config,
    const pair<size_t, size_t> &lattice_id_jump_pair)
{
  // Get Pair Direction
  Vector3d
      relative_distance_vector_pair =
          config.GetRelativeDistanceVectorLattice(lattice_id_jump_pair.first,
                                                  lattice_id_jump_pair.second);

  const Matrix3d basis = config.GetBasis();

  RowVector3d
      pair_direction = relative_distance_vector_pair.transpose() * basis;
  pair_direction.normalize();

  RowVector3d vertical_vector;

  // First nearest neighbors
  vector<size_t> nn_list =
      config.GetNeighborLatticeIdVectorOfLattice(lattice_id_jump_pair.first, 1);

  sort(nn_list.begin(), nn_list.end());

  for (const auto nn_id : nn_list)
  {
    // Get Jump Vector
    Vector3d
        relative_distance_vector =
            config.GetRelativeDistanceVectorLattice(lattice_id_jump_pair.first, nn_id);

    RowVector3d jump_vector = relative_distance_vector.transpose() * basis;
    jump_vector.normalize();

    // Check for vertical direction
    const double dot_prod = pair_direction.dot(jump_vector);
    if (abs(dot_prod) < constants::kEpsilon)
    {
      vertical_vector = jump_vector;
      break;
    }
  }

  // Rotation Matrix
  Matrix3d rotation_matrix;
  rotation_matrix.row(0) = pair_direction;
  rotation_matrix.row(1) = vertical_vector;
  rotation_matrix.row(2) = pair_direction.cross(vertical_vector);

  return rotation_matrix.transpose();
}

vector<vector<size_t>> GetEquivalentSites3Fold(const Config &config,
                                               const size_t maxBondOrder)
{
  // Does not depends on jump pairs.
  pair<size_t, size_t> jumpPair = {0, config.GetNeighborLatticeIdVectorOfLattice(0, 1)[0]};

  auto symmetricallySortedVector = GetSortedLatticeVectorStateOfPair(config,
                                                                     jumpPair,
                                                                     maxBondOrder);

  vector<size_t> equivalentSites;
  vector<vector<size_t>> equivalentSiteVector;

  size_t idx = 0;

  while (idx < symmetricallySortedVector.size())
  {
    auto latticeId = symmetricallySortedVector[idx];

    auto directionVectorI =
        config.GetNormalizedDirection(latticeId, jumpPair.first);
    auto directionVectorII =
        config.GetNormalizedDirection(latticeId, jumpPair.second);

    if (directionVectorI == directionVectorII)
    {
      equivalentSites.emplace_back(idx);
      idx++;
    }
    else
    {
      // Sample in threes
      // Works only till 2nd NN
      for (size_t id = idx; id < idx + 3; id++)
      {
        equivalentSites.emplace_back(id);
      }
      idx += 3;
    }

    equivalentSiteVector.emplace_back(equivalentSites);
    equivalentSites = {};
  }

  return equivalentSiteVector;
}

/*
// Pairs between Migrating Atom and Nearest Neigbours
VectorXd GetEncodingMigratingAtomPair (
  const Config &config,
  const vector<vector<size_t>> &equivalentSitesEncoding,
  const vector<size_t> &symmetricallySortedVector,
  const unordered_map<string, RowVectorXd> &oneHotEncodingMap,
  const Element &migratingAtom)
{

// maxBondOrder is used for getting the NN for the jumpPair
// Currently only support till maxBondOrder = 2

// for bondOrder 1, there are 6 equivalent groups
// the two central groups will be the first nn of the migrating atom (6 atoms)
// then the next to each of the central groups will be the second nn atom
// for the migrating atom at transition site
// and third nn for the migrating atom based on the cutoff described
// so in the first nn shell of the jump pair, there will be 14 nn of the migrating Atom

  // Pair of Migrating Atom and the NN
  // Not required storing the element pairs for verification
  // vector<vector<string>> equivalentPairsClusters.reserve();

  size_t sizeEncodeNonSymmPairs = 4;

  VectorXd encodeVector;

  for (auto &equivalentSites : equivalentSitesEncoding)
  {
    vector<string> equivalentPairs;

    // Since non symmetric pairs
    // Migrating Atom will be first site
    // Need to work on this so as to generalize this

    RowVectorXd pairEncodeVector = RowVectorXd::Zero(Index(sizeEncodeNonSymmPairs));

    for (auto &sites : equivalentSites)
    {
      auto latticeId = symmetricallySortedVector[sites];
      auto neighborElement = config.GetElementOfLattice(latticeId);

      string elementPair = migratingAtom.GetElementString() +
                                neighborElement.GetElementString();

      equivalentPairs.emplace_back(elementPair);

      RowVectorXd oneHotVector = oneHotEncodingMap.at(elementPair);

      pairEncodeVector += oneHotVector;
    }

  // equivalentPairsClusters.emplace_back(equivalentPairs);

  encodeVector.conservativeResize(encodeVector.size() + pairEncodeVector.size());

  // Concatenate pairEncodeVector to the resized vector
  encodeVector.tail(pairEncodeVector.size()) = pairEncodeVector;

  }

  // print2DStringVector(equivalentPairsClusters);

  return encodeVector;
}

*/

VectorXd GetEncodingMigratingAtomPair(
    const Config &config,
    const vector<vector<size_t>> &equivalentSitesEncoding,
    const vector<size_t> &symmetricallySortedVector,
    const unordered_map<string, RowVectorXd> &oneHotEncodingMap,
    const Element &migratingAtom)
{
  // Size of non-symmetric pairs
  size_t numElements = 2; // Considering only binary for now
  size_t sizeEncodeNonSymmPairs = pow(numElements, 2);

  // Calculate the total size needed for encodeVector
  size_t totalSize;
  totalSize = equivalentSitesEncoding.size() * sizeEncodeNonSymmPairs;

  // Pre-allocate encodeVector with the total expected size
  VectorXd encodeVector(totalSize);
  size_t offset = 0;
  
  // cout << totalSize << endl;
  // print2DVector(equivalentSitesEncoding);
  
  // Loop through the equivalent sites encoding
  for (const auto &equivalentSites : equivalentSitesEncoding)
  {
    RowVectorXd pairEncodeVector = RowVectorXd::Zero(sizeEncodeNonSymmPairs);
    
    // Loop through the equivalent sites and build the pair encoding vector
    for (auto &sites : equivalentSites)
    {
      auto latticeId = symmetricallySortedVector[sites];
      auto neighborElement = config.GetElementOfLattice(latticeId);

      string elementPair = migratingAtom.GetElementString() +
                           neighborElement.GetElementString();

      // cout << elementPair << endl;

      try
      {
        pairEncodeVector += oneHotEncodingMap.at(elementPair); 
      }
      catch (const std::out_of_range &e)
      {
        for (auto latticeId : symmetricallySortedVector)
        {
          cout << latticeId << config.GetElementOfLattice(latticeId) << " ";
        }
        cout << endl;
        std::cout << "Error: Missing Element Pair for " << elementPair << "(" << latticeId << ")" << std::endl;
        exit(1);
      }
    }

    // Store the result into encodeVector at the correct offset
    encodeVector.segment(offset, pairEncodeVector.size()) = pairEncodeVector;
    offset += pairEncodeVector.size(); // Move the offset for next pair
  }

  // Return the final accumulated encode vector
  return encodeVector;
}
