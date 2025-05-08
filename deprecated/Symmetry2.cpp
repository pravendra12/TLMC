#include "Symmetry.h"

/**
 * @brief Finds the closest matching index for a target vector in a map of positions to indices.
 *
 * Iterates through a map of 3D vectors (positions) and indices, calculating the Euclidean distance
 * between each position and the target vector. If the distance is within a predefined tolerance
 * (constants::kEpsilon), the corresponding index is returned. If no match is found, an invalid
 * index (numeric_limits<size_t>::max()) is returned.
 *
 * @param positionToIndex Map of 3D vectors (positions) to indices.
 * @param target Target 3D vector.
 * @return Closest matching index or numeric_limits<size_t>::max().
 */
static size_t findClosestMatch(
    const unordered_map<Vector3d, size_t, Vector3dHash> &positionToIndex,
    const Vector3d &target)
{

  for (const auto &[pos, index] : positionToIndex)
  {
    if ((pos - target).norm() < constants::kEpsilon) // Check if within tolerance
    {
      return index; // Return the matching index
    }
  }
  return numeric_limits<size_t>::max(); // Return an invalid index if no match is found
}

// Convert degrees to radians
inline double toRadians(double degrees)
{
  return degrees * M_PI / 180.0;
}

/**
 * @brief Rotates a 3D point around a specified axis by a given angle, with
 * respect to a center point. This function uses Rodrigues' rotation formula
 * to perform the rotation.
 *
 * @param point The 3D point to be rotated (Vector3d).
 * @param axis The axis of rotation, represented as a 3D vector (Vector3d).
 * It will be normalized internally.
 * @param theta The angle of rotation in degrees.
 * @param center The center point for rotation (Vector3d).
 * @return The rotated 3D point (Vector3d).
 */
static Vector3d rotatePoints(const Vector3d &point,
                             const Vector3d &axis,
                             double theta,
                             const Vector3d &center)
{

  // Normalize axis
  Vector3d k = axis.normalized();
  double theta_rad = toRadians(theta);
  double cos_theta = cos(theta_rad);
  double sin_theta = sin(theta_rad);

  // Translate point relative to center
  Vector3d p = point - center;

  // Rodrigues' rotation formula: p_rot = p cosθ + (k × p) sinθ + k (k · p) (1 - cosθ)
  Vector3d p_rot = p * cos_theta +
                   k.cross(p) * sin_theta +
                   k * (k.dot(p)) * (1.0 - cos_theta);

  // Translate back
  return p_rot + center;
}

/**
 * @brief Computes the center position of a lattice pair in a periodic lattice system.
 *
 * Calculates the center position of two lattice points, considering periodic boundary
 * conditions. Adjusts the relative positions of the lattice points based on their
 * periodicity to ensure the correct center is computed.
 *
 * @param config Configuration object providing access to lattice positions.
 * @param lattice_id_jump_pair Pair of lattice IDs representing the two lattice points.
 * @return RowVector3d representing the center position of the lattice pair.
 */
static RowVector3d GetLatticePairCenter(
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

/**
 * @brief Finds the closest <111> direction to a given vector.
 *
 * This function determines the closest <111> crystallographic direction
 * to the input vector by calculating the dot product with all possible
 * <111> directions and selecting the one with the maximum absolute value.
 * If multiple directions have the same dot product, a tie-breaking rule
 * based on component comparison is applied.
 *
 * @param pairDirection The input vector to compare against <111> directions.
 * @return Vector3d The closest <111> direction as a normalized vector.
 */
static Vector3d GetClosest111Direction(const Vector3d &pairDirection)
{
  Vector3d normDir = pairDirection.normalized();

  vector<Vector3d> directions = {
      Vector3d(1, 1, 1).normalized(),
      Vector3d(-1, 1, 1).normalized(),
      Vector3d(1, -1, 1).normalized(),
      Vector3d(-1, -1, 1).normalized(),
      Vector3d(1, 1, -1).normalized(),
      Vector3d(-1, 1, -1).normalized(),
      Vector3d(1, -1, -1).normalized(),
      Vector3d(-1, -1, -1).normalized()};

  Vector3d closestDir = directions[0];

  // Return the absolute value
  double maxDot = fabs(normDir.dot(closestDir));

  for (const auto &dir : directions)
  {
    double dot = fabs(normDir.dot(dir));
    if (dot > maxDot + constants::kEpsilon)
    {
      maxDot = dot;
      closestDir = dir;
    }
    else if (fabs(dot - maxDot) < constants::kEpsilon)
    {
      for (int i = 0; i < 3; ++i)
      {
        if (fabs(dir[i]) > constants::kEpsilon &&
            fabs(closestDir[i]) > constants::kEpsilon)
        {
          if (dir[i] > closestDir[i])
          {
            closestDir = dir;
            maxDot = dot;
          }
          break;
        }
      }
    }
  }

  return closestDir;
}

/**
 * @brief Generates a set of points equivalent to a given point under 3-bar symmetry.
 *
 * This function computes the equivalent points of a given starting point under
 * 3-bar symmetry (a combination of 3-fold rotational symmetry and inversion symmetry).
 * The symmetry is applied around a specified axis and center.
 *
 * @param start_point The initial point to generate equivalent points for.
 * @param axis The axis of symmetry for the 3-bar symmetry operation.
 * @param center The center of symmetry for the 3-bar symmetry operation.
 * @return A vector of Vector3d objects representing the equivalent points.
 *
 * The function performs the following operations:
 * 1. Rotates the starting point around the axis by angles of 0°, 120°, and 240°,
 *    adding the resulting points to the equivalent points list.
 * 2. Applies inversion symmetry to the rotated points (relative to the center)
 *    and adds these inverted points to the equivalent points list.
 *
 * The resulting vector contains six points: three from rotation and three from
 * inversion symmetry.
 */
static vector<Vector3d> GetEquivalentPositionUnder3BarSymmetry(
    const Vector3d &start_point,
    const Vector3d &axis,
    const Vector3d &center)
{
  vector<Vector3d> equivalent_points;
  equivalent_points.reserve(6);

  vector<double> angles = {0.0, 120.0, 240.0};

  for (double angle : angles)
  {
    equivalent_points.push_back(rotatePoints(start_point, axis, angle, center));
  }

  for (double angle : angles)
  {
    Vector3d rotated = rotatePoints(start_point, axis, angle, center);
    equivalent_points.push_back(2.0 * center - rotated);
  }

  return equivalent_points;
}

vector<size_t> GetSortedLatticeStatesForPairUnder3BarSymmetry(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair,
    const size_t &maxBondOrder)
{
  auto neighboring_lattice_ids = config.GetNeighboringLatticeIdSetOfPair(latticeIdJumpPair,
                                                                         maxBondOrder);
  size_t num_sites = neighboring_lattice_ids.size();

  Vector3d first_pos = config.GetCartesianPositionOfLattice(latticeIdJumpPair.first);
  Vector3d second_pos = config.GetCartesianPositionOfLattice(latticeIdJumpPair.second);
  Vector3d pair_direction = first_pos - second_pos;
  Vector3d rotation_axis = GetClosest111Direction(pair_direction);

  // Expects cartesian position
  RowVector3d center = GetLatticePairCenter(config, latticeIdJumpPair);
  center = center * config.GetBasis();

  Vector3d z_axis(0, 0, 1);
  Matrix3d rotation_matrix;
  if (rotation_axis.dot(z_axis) < 1.0 - constants::kEpsilon)
  {
    Vector3d cross = rotation_axis.cross(z_axis);
    double angle = acos(rotation_axis.dot(z_axis));
    rotation_matrix = Eigen::AngleAxisd(angle, cross.normalized()).toRotationMatrix();
  }
  else
  {
    rotation_matrix.setIdentity();
  }

  vector<tuple<size_t, double, RowVector3d>> lattice_id_keys;
  lattice_id_keys.reserve(num_sites);
  for (const auto id : neighboring_lattice_ids)
  {
    RowVector3d pos = config.GetCartesianPositionOfLattice(id).transpose();
    pos -= center;
    pos = pos * rotation_matrix;
    double key = pos.norm();
    lattice_id_keys.emplace_back(id, key, pos);
  }

  sort(lattice_id_keys.begin(), lattice_id_keys.end(), [](const auto &lhs, const auto &rhs)
       {
          double key_diff = get<1>(lhs) - get<1>(rhs);
          if (fabs(key_diff) > constants::kEpsilon) 
          {
            return key_diff < 0;
          }
          const auto& pos_lhs = get<2>(lhs);
          const auto& pos_rhs = get<2>(rhs);
          for (int i = 0; i < 3; ++i) {
          double diff = pos_lhs[i] - pos_rhs[i];
          if (abs(diff) > constants::kEpsilon) 
          {
            return diff < 0;
          }
          }
          return false; });

  vector<size_t> sorted_lattice_ids;
  sorted_lattice_ids.reserve(num_sites);
  for (const auto &item : lattice_id_keys)
  {
    sorted_lattice_ids.push_back(get<0>(item));
  }
  return sorted_lattice_ids;
}

vector<vector<size_t>> GetEquivalentSitesUnder3BarSymmetry(
    const Config &config,
    size_t maxBondOrder,
    const pair<size_t, size_t> &latticeIdPair)
{
  auto ssVector = GetSortedLatticeStatesForPairUnder3BarSymmetry(config,
                                                        latticeIdPair,
                                                        maxBondOrder);

  Vector3d centralPos = config.GetCartesianPositionOfLattice(latticeIdPair.first);
  Vector3d nnPos = config.GetCartesianPositionOfLattice(latticeIdPair.second);
  Vector3d pairDirection = centralPos - nnPos;
  Vector3d rotationAxis = GetClosest111Direction(pairDirection);

  Vector3d center = GetLatticePairCenter(config, latticeIdPair);

  unordered_map<Vector3d, size_t, Vector3dHash> positionToIndex;
  vector<Vector3d> cartesianPositionVector;
  cartesianPositionVector.reserve(ssVector.size());

  for (size_t i = 0; i < ssVector.size(); ++i)
  {
    auto pos = config.GetCartesianPositionOfLattice(ssVector[i]);
    cartesianPositionVector.push_back(pos);
    positionToIndex[pos] = i;
  }

  vector<bool> processed(ssVector.size(), false);
  vector<vector<size_t>> equivalentEncodingVector;
  for (size_t i = 0; i < ssVector.size(); ++i)
  {
    if (processed[i])
      continue;

    vector<size_t> equivalentIndices = {i};
    processed[i] = true;

    // Get equivalent points with inversion
    /*
    double theta = 60.0; // 3 Bar Symmetry
    auto equivPositions = GetEquivalentPoints(cartesianPositionVector[i],
                                              rotationAxis,
                                              theta,
                                              center,
                                              true);

    */
    auto equivPositions = GetEquivalentPositionUnder3BarSymmetry(cartesianPositionVector[i],
                                                                 rotationAxis,
                                                                 center);

    for (const auto &equivPos : equivPositions)
    {
      size_t index = findClosestMatch(positionToIndex, equivPos);
      if (index != numeric_limits<size_t>::max() && !processed[index])
      {
        equivalentIndices.push_back(index);
        processed[index] = true;
      }
    }

    sort(equivalentIndices.begin(), equivalentIndices.end());
    equivalentEncodingVector.push_back(move(equivalentIndices));
  }

  sort(equivalentEncodingVector.begin(), equivalentEncodingVector.end(),
       [](const auto &a, const auto &b)
       { return a[0] < b[0]; });

  // cout  << "Equivalent Sites for Pair " << latticeIdPair.first << " -> " << latticeIdPair.second << ": ";
  //  print2DVector(equivalentEncodingVector);

  return equivalentEncodingVector;
}

/*
vector<Vector3d> GetEquivalentPoints(const Vector3d &start_point,
                                     const Vector3d &axis,
                                     double theta,
                                     const Vector3d &center)
{
  int n = 360.0 / theta;

  // cout << "number of rotations: " << n << endl;

  vector<Vector3d> equivalent_points;
  equivalent_points.reserve(n); // Pre-allocate space

  for (int i = 0; i < n; ++i)
  {
    double angle = i * theta;
    Vector3d rotated = rotatePoints(start_point, axis, angle, center);
    equivalent_points.push_back(rotated);
    // cout << rotated.transpose() << endl;
  }

  return equivalent_points;
}
  */

/*
vector<Vector3d> GetEquivalentPoints(const Vector3d &startPoint,
                                     const Vector3d &axis,
                                     double theta,
                                     const Vector3d &center,
                                     bool applyInversion)
{
  int n = static_cast<int>(std::round(360.0 / theta));

  vector<Vector3d> equivalentPoints;
  equivalentPoints.reserve(n);

  for (int i = 0; i < n; ++i)
  {
    double angle = i * theta;
    Vector3d rotated = rotatePoints(startPoint, axis, angle, center);

    if (applyInversion)
    {
      rotated = 2.0 * center - rotated; // Inversion through the center
    }

    equivalentPoints.push_back(rotated);
  }

  return equivalentPoints;
}
*/

/*
vector<vector<size_t>> GetEquivalentSitesUnderKBarSymmetry(const Config &config,
                                                           size_t maxBondOrder,
                                                           size_t kFoldRotation)
{
  vector<vector<size_t>> equivalentEncodingVectorKBar;

  // Get lattice pair (central and nearest neighbor)
  const size_t centralLatticeId = config.GetCentralAtomLatticeId();
  const auto nnLatticeIdVector = config.GetNeighborLatticeIdVectorOfLattice(centralLatticeId, 1);

  const size_t nnLatticeId = nnLatticeIdVector[0];
  const pair<size_t, size_t> latticeIdPair = {centralLatticeId, nnLatticeId};

  // Get symmetrically sorted lattice ID vector for the pair
  const auto ssVector = config.GetSortedLatticeVectorStateOfPair(latticeIdPair, maxBondOrder);

  // Transition position
  const Vector3d centralLatticePosition = config.GetCartesianPositionOfLattice(centralLatticeId);
  const Vector3d nnLatticePosition = config.GetCartesianPositionOfLattice(nnLatticeId);
  const Vector3d transitionPosition = 0.5 * (centralLatticePosition + nnLatticePosition);

  // Rotation axis
  const Vector3d rotationAxis = centralLatticePosition - nnLatticePosition;
  // cout << "Transition Position: " << transitionPosition.transpose() << endl;

  // Pre-populate position vector
  vector<Vector3d> cartesianPositionVector;
  cartesianPositionVector.reserve(ssVector.size());

  for (const auto &latticeId : ssVector)
  {
    cartesianPositionVector.emplace_back(config.GetCartesianPositionOfLattice(latticeId));
  }

  // Use a hash map for O(1) lookup of positions to indices
  unordered_map<Vector3d, size_t, Vector3dHash> positionToIndex;
  for (size_t i = 0; i < cartesianPositionVector.size(); ++i)
  {
    positionToIndex[cartesianPositionVector[i]] = i;
    // cout << cartesianPositionVector[i].transpose() << " : " << i << endl;
  }

  auto equivalentEncodingVectorKFold = GetEquivalentSitesUnderKFoldRotation(config,
                                                                            maxBondOrder,
                                                                            kFoldRotation);
  unordered_map<size_t, size_t> invertedSiteMap;

  for (size_t i = 0; i < equivalentEncodingVectorKFold.size(); ++i)
  {

    // set<size_t> allIndices(equivalentEncodingVector[i].begin(), equivalentEncodingVector[i].end());

    // Invert each position in group i
    for (auto id : equivalentEncodingVectorKFold[i])
    {
      Vector3d pos = config.GetCartesianPositionOfLattice(ssVector[id]); // individual site position
      Vector3d inverted = 2.0 * transitionPosition - pos;

      // Find the closest matching representative position
      size_t matchedId = findClosestMatch(positionToIndex, inverted);

      // cout << id << " " << pos.transpose() << " Inverted " << inverted.transpose() << " " << matchedId << endl;

      invertedSiteMap.insert(make_pair(id, matchedId));
    }
  }

  // combine the inversion and 3 bar encodings

  std::unordered_set<size_t> visited;

  for (const auto &equivalentSites : equivalentEncodingVectorKFold)
  {
    std::vector<size_t> combinedEquivalentSites;
    for (auto siteId : equivalentSites)
    {
      if (visited.count(siteId))
        continue;

      combinedEquivalentSites.push_back(siteId);
      visited.insert(siteId);

      if (invertedSiteMap.count(siteId))
      {
        size_t pairId = invertedSiteMap[siteId];
        if (!visited.count(pairId))
        {
          combinedEquivalentSites.push_back(pairId);
          visited.insert(pairId);
        }
      }
    }

    if (!combinedEquivalentSites.empty())
      equivalentEncodingVectorKBar.push_back(combinedEquivalentSites);
  }

  // print2DVector(equivalentEncodingVectorKBar);

  return equivalentEncodingVectorKBar;
}

vector<vector<size_t>> GetEquivalentSitesUnderKFoldRotation(const Config &config,
                                                            size_t maxBondOrder,
                                                            size_t kFoldRotation)
{
  // Get lattice pair (central and nearest neighbor)
  // const size_t centralLatticeId = config.GetCentralAtomLatticeId();
  const size_t centralLatticeId = config.GetCentralAtomLatticeId();

  const auto nnLatticeIdVector = config.GetNeighborLatticeIdVectorOfLattice(centralLatticeId, 1);

  const size_t nnLatticeId = nnLatticeIdVector[0];
  const pair<size_t, size_t> latticeIdPair = {centralLatticeId, nnLatticeId};

  // Get symmetrically sorted lattice ID vector for the pair
  const auto ssVector = config.GetSortedLatticeVectorStateOfPair(latticeIdPair, maxBondOrder);

  // Transition position
  const Vector3d centralLatticePosition = config.GetCartesianPositionOfLattice(centralLatticeId);
  const Vector3d nnLatticePosition = config.GetCartesianPositionOfLattice(nnLatticeId);
  const Vector3d transitionPosition = 0.5 * (centralLatticePosition + nnLatticePosition);

  // Rotation axis
  const Vector3d rotationAxis = centralLatticePosition - nnLatticePosition;
  // cout << "Transition Position: " << transitionPosition.transpose() << endl;

  // Pre-populate position vector
  vector<Vector3d> cartesianPositionVector;
  cartesianPositionVector.reserve(ssVector.size());

  for (const auto &latticeId : ssVector)
  {
    cartesianPositionVector.emplace_back(config.GetCartesianPositionOfLattice(latticeId));
  }

  // Store equivalent sites and their encodings
  vector<vector<size_t>> equivalentEncodingVector;
  // vector<Vector3d> equivalentEncodingPositionVector;

  // Use a hash map for O(1) lookup of positions to indices
  unordered_map<Vector3d, size_t, Vector3dHash> positionToIndex;
  for (size_t i = 0; i < cartesianPositionVector.size(); ++i)
  {
    positionToIndex[cartesianPositionVector[i]] = i;
    // cout << cartesianPositionVector[i].transpose() << " : " << i << endl;
  }

  double rotationAngle = 360.0 / kFoldRotation; // 120 degrees for 3-fold

  vector<bool> processed(cartesianPositionVector.size(), false);

  for (size_t i = 0; i < cartesianPositionVector.size(); ++i)
  {
    if (processed[i])
      continue;

    const Vector3d &position = cartesianPositionVector[i];
    auto equivalentPositions = GetEquivalentPoints(position, rotationAxis, rotationAngle, transitionPosition, false);

    vector<size_t> equivalentIndices = {i}; // Include the original position
    processed[i] = true;

    // Check equivalent positions

    for (const auto &equivPos : equivalentPositions)
    {
      size_t index = findClosestMatch(positionToIndex, equivPos);

      if (index != numeric_limits<size_t>::max()) // If a valid index is found
      {
        if (!processed[index])
        {
          equivalentIndices.push_back(index);
          processed[index] = true;
        }
      }
    }

    equivalentEncodingVector.push_back(move(equivalentIndices));
    // equivalentEncodingPositionVector.push_back(position);
  }

  // roto inversion

  for (size_t i = 0; i < equivalentEncodingVector.size(); ++i)
  {

    // set<size_t> allIndices(equivalentEncodingVector[i].begin(), equivalentEncodingVector[i].end());

    // Invert each position in group i
    for (auto id : equivalentEncodingVector[i])
    {
      Vector3d pos = config.GetCartesianPositionOfLattice(ssVector[id]); // individual site position
      Vector3d inverted = 2.0 * transitionPosition - pos;

      // Find the closest matching representative position
      size_t matchedId = findClosestMatch(positionToIndex, inverted);

      cout << id << " " << pos.transpose() << " Inverted " << inverted.transpose() << " " << matchedId << endl;
    }
  }

  // print2DVector(equivalentEncodingVector);

  return equivalentEncodingVector;
}
*/

/*
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
*/

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
/*
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

      try
      {
        pairEncodeVector += oneHotEncodingMap.at(elementPair);
      }
      catch (const out_of_range &e)
      {
        cout << "Error: Missing Element Pair for " << elementPair << "(" << latticeId << ")" << endl;
        exit(1);
      }
    }
    // Normalized encoding vector
    // Store the result into encodeVector at the correct offset
    encodeVector.segment(offset, pairEncodeVector.size()) = pairEncodeVector / equivalentSites.size();
    offset += pairEncodeVector.size(); // Move the offset for next pair
  }

  // Return the final accumulated encode vector
  return encodeVector;
}
*/
