#include "Symmetry.h"

const double K_EPSILON = 0.01;

/**
 * @brief Finds the closest matching index for a target vector in a map of positions to indices.
 *
 * Iterates through a map of 3D vectors (positions) and indices, calculating the Euclidean distance
 * between each position and the target vector. If the distance is within a predefined tolerance
 * (K_EPSILON ), the corresponding index is returned. If no match is found, an invalid
 * index (numeric_limits<size_t>::max()) is returned.
 *
 * @param positionToIndex Map of 3D vectors (positions) to indices.
 * @param target Target 3D vector.
 * @return Closest matching index or numeric_limits<size_t>::max().

static size_t findClosestMatch(
    const unordered_map<Vector3d, size_t, Vector3dHash> &positionToIndex,
    const Vector3d &target)
{

  for (const auto &[pos, index] : positionToIndex)
  {
    if ((pos - target).norm() < K_EPSILON ) // Check if within tolerance
    {
      return index; // Return the matching index
    }
  }
  return numeric_limits<size_t>::max(); // Return an invalid index if no match is found
}
  */

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
/*
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
    if (dot > maxDot + K_EPSILON )
    {
      maxDot = dot;
      closestDir = dir;
    }
    else if (fabs(dot - maxDot) < K_EPSILON )
    {
      for (int i = 0; i < 3; ++i)
      {
        if (fabs(dir[i]) > K_EPSILON  &&
            fabs(closestDir[i]) > K_EPSILON )
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
  */

Vector3d GetClosest111Direction(const Vector3d &direction)
{
  static const vector<Vector3d> k111Directions = {
      {1, 1, 1}, {-1, 1, 1}, {1, -1, 1}, {1, 1, -1}, {-1, -1, 1}, {-1, 1, -1}, {1, -1, -1}, {-1, -1, -1}};

  Vector3d norm_direction = direction.normalized();
  Vector3d best_direction = k111Directions[0];
  double max_dot = -1.0;

  for (const auto &dir : k111Directions)
  {
    double dot = fabs(norm_direction.dot(dir.normalized()));
    if (dot > max_dot + K_EPSILON)
    {
      max_dot = dot;
      best_direction = dir;
    }
    else if (fabs(dot - max_dot) < K_EPSILON)
    {
      // Break ties lexicographically
      if (dir(0) < best_direction(0) ||
          (dir(0) == best_direction(0) && dir(1) < best_direction(1)) ||
          (dir(0) == best_direction(0) && dir(1) == best_direction(1) && dir(2) < best_direction(2)))
      {
        best_direction = dir;
      }
    }
  }

  return best_direction.normalized();
}

/**
 * @brief Computes the relative position of a lattice pair center in a periodic lattice system.
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
 * @brief Wraps a 3D vector into the periodic unit cell [0, 1).
 *
 * @param pos The 3D vector to wrap.
 * @return Vector3d The wrapped position.
 */
#include <cmath>

static Vector3d wrapToUnitCell(const Vector3d &pos)
{
  Vector3d wrapped;
  for (int i = 0; i < 3; ++i)
  {
    wrapped[i] = pos[i] - std::floor(pos[i]);
  }
  return wrapped;
}

/**
 * @brief Computes the periodic distance between two 3D vectors.
 *
 * Accounts for periodic boundary conditions by considering the minimum distance
 * across periodic images.
 *
 * @param a First 3D vector.
 * @param b Second 3D vector.
 * @return double The minimum Euclidean distance.
 */
static double periodicDistance(const Vector3d &a, const Vector3d &b)
{
  Vector3d diff = a - b;
  for (int i = 0; i < 3; ++i)
  {
    diff[i] = diff[i] - round(diff[i]); // Find the shortest periodic image
  }
  return diff.norm();
}

/**
 * @brief Finds the closest matching index for a target vector in a map of positions to indices.
 *
 * Iterates through a map of 3D vectors (positions) and indices, calculating the periodic
 * Euclidean distance between each position and the target vector. If the distance is
 * within a predefined tolerance (K_EPSILON ), the corresponding index is returned.
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
    if (periodicDistance(pos, target) < K_EPSILON)
    {
      return index;
    }
  }
  return numeric_limits<size_t>::max();
}

/**
 * @brief Generates a set of points equivalent to a given point under 3-bar symmetry.
 *
 * Computes equivalent points under 3-fold rotational symmetry and inversion symmetry,
 * ensuring all positions are wrapped into the unit cell.
 *
 * @param start_point The initial point to generate equivalent points for.
 * @param axis The axis of symmetry for the 3-bar symmetry operation.
 * @param center The center of symmetry for the 3-bar symmetry operation.
 * @return A vector of Vector3d objects representing the equivalent points.
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
    Vector3d rotated = rotatePoints(start_point, axis, angle, center);
    equivalent_points.push_back(wrapToUnitCell(rotated));
  }

  for (double angle : angles)
  {
    Vector3d rotated = rotatePoints(start_point, axis, angle, center);
    Vector3d inverted = 2.0 * center - rotated;
    equivalent_points.push_back(wrapToUnitCell(inverted));
  }

  cout << "Equivalent Points: " << endl;
  for (auto vec : equivalent_points)
  {
    cout << vec.transpose() << endl;
  }

  return equivalent_points;
}

/**
 * @brief Groups lattice sites into equivalent sets under 3-bar symmetry.
 *
 * Uses the sorted lattice IDs and computes equivalent positions, handling periodic
 * boundary conditions to ensure correct grouping for boundary pairs.
 *
 * @param config Configuration object providing lattice positions.
 * @param maxBondOrder Maximum bond order for neighbor search.
 * @param latticeIdPair Pair of lattice IDs.
 * @return Vector of vectors containing indices of equivalent sites.
 */
vector<vector<size_t>> GetEquivalentSitesUnder3BarSymmetry(
    const Config &config,
    size_t maxBondOrder,
    const pair<size_t, size_t> &latticeIdPair)
{
  // Get sorted lattice IDs
  auto ssVector = GetSortedLatticeStatesForPairUnder3BarSymmetry(config, latticeIdPair, maxBondOrder);

  // Compute pair direction and rotation axis
  Vector3d centralPos = config.GetRelativePositionOfLattice(latticeIdPair.first);
  Vector3d nnPos = config.GetRelativePositionOfLattice(latticeIdPair.second);
  Vector3d pairDirection = centralPos - nnPos;
  Vector3d rotationAxis = GetClosest111Direction(pairDirection);

  // Compute center
  Vector3d center = GetLatticePairCenter(config, latticeIdPair).transpose();

  // Store positions and indices
  unordered_map<Vector3d, size_t, Vector3dHash> positionToIndex;
  vector<Vector3d> relativePositionVector;
  relativePositionVector.reserve(ssVector.size());

  for (size_t i = 0; i < ssVector.size(); ++i)
  {
    auto pos = config.GetRelativePositionOfLattice(ssVector[i]);
    relativePositionVector.push_back(pos);
    positionToIndex[pos] = i;
  }

  // Track processed sites
  vector<bool> processed(ssVector.size(), false);
  vector<vector<size_t>> equivalentEncodingVector;

  for (size_t i = 0; i < ssVector.size(); ++i)
  {
    if (processed[i])
    {
      continue;
    }

    vector<size_t> equivalentIndices = {i};
    processed[i] = true;

    // Compute equivalent positions with PBC
    auto equivPositions = GetEquivalentPositionUnder3BarSymmetry(
        relativePositionVector[i], rotationAxis, center);

    for (const auto &equivPos : equivPositions)
    {
      size_t index = findClosestMatch(positionToIndex, equivPos);
      if (index != numeric_limits<size_t>::max() && !processed[index])
      {
        equivalentIndices.push_back(index);
        processed[index] = true;
      }
    }

    // Sort indices within each equivalent set
    sort(equivalentIndices.begin(), equivalentIndices.end());
    equivalentEncodingVector.push_back(move(equivalentIndices));
  }

  // Sort equivalent sets by the first index
  sort(equivalentEncodingVector.begin(), equivalentEncodingVector.end(),
       [](const auto &a, const auto &b)
       { return a[0] < b[0]; });

  return equivalentEncodingVector;
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

inline Eigen::RowVector3d C3(const Eigen::RowVector3d &v)
{
  return Eigen::RowVector3d(-v[1], v[0] - v[1], v[2]); // 120° rotation
}

inline Eigen::RowVector3d C3Sq(const Eigen::RowVector3d &v)
{
  return C3(C3(v)); // 240°
}

inline Eigen::RowVector3d Mod1(const Eigen::RowVector3d &v)
{
  return (v.array() - v.array().floor()).matrix(); // wrap to [0,1)
}

inline Eigen::RowVector3d CanonicalUnder3Bar(const Eigen::RowVector3d &v)
{
  std::array<Eigen::RowVector3d, 6> images = {
      Mod1(v),
      Mod1(C3(v)),
      Mod1(C3Sq(v)),
      Mod1(-v),
      Mod1(-C3(v)),
      Mod1(-C3Sq(v))};

  return *std::min_element(images.begin(), images.end(), [](const auto &a, const auto &b)
                           {
    for (int i = 0; i < 3; ++i) {
      if (a[i] + K_EPSILON  < b[i]) return true;
      if (a[i] - K_EPSILON  > b[i]) return false;
    }
    return false; });
}

inline bool PositionCompare3BarSymmetry(
    const std::pair<size_t, Eigen::RowVector3d> &lhs,
    const std::pair<size_t, Eigen::RowVector3d> &rhs)
{

  const auto lhs_canon = CanonicalUnder3Bar(lhs.second);
  const auto rhs_canon = CanonicalUnder3Bar(rhs.second);

  for (int i = 0; i < 3; ++i)
  {
    const double diff = lhs_canon[i] - rhs_canon[i];
    if (diff < -K_EPSILON)
      return true;
    if (diff > K_EPSILON)
      return false;
  }

  return false;
}

vector<size_t> GetSortedLatticeStatesForPairUnder3BarSymmetry(const Config &config,
                                                              const pair<size_t, size_t> &lattice_id_jump_pair,
                                                              const size_t &max_bond_order)
{
  auto neighboring_lattice_ids = config.GetNeighboringLatticeIdSetOfPair(lattice_id_jump_pair,
                                                                         max_bond_order);

  size_t num_sites = neighboring_lattice_ids.size();
  Eigen::RowVector3d move_distance =
      Eigen::RowVector3d(0.5, 0.5, 0.5) - GetLatticePairCenter(config, lattice_id_jump_pair);

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

  // Rotate lattice vectors
  // Config::RotateLatticeVector(lattice_id_hashmap, config.GetLatticePairRotationMatrix(lattice_id_jump_pair));

  // Convert unordered_map to vector for sorting
  std::vector<std::pair<size_t, Eigen::RowVector3d>>
      lattice_id_vector(lattice_id_hashmap.begin(), lattice_id_hashmap.end());

  // Sort the lattice vector based on PositionCompare
  std::sort(lattice_id_vector.begin(), lattice_id_vector.end(), PositionCompare3BarSymmetry);

  // Extract and return only the lattice IDs
  std::vector<size_t> sorted_lattice_ids;
  sorted_lattice_ids.reserve(lattice_id_vector.size());
  for (const auto &pair : lattice_id_vector)
  {
    sorted_lattice_ids.push_back(pair.first);
    // std::cout << pair.first << " " << pair.second << std::endl;
  }
  // cout  << "Jump Pair: " << lattice_id_jump_pair.first << " -> " << lattice_id_jump_pair.second << endl;
  // cout  << "Rotation Axis: " << rotation_axis.transpose() << endl;
  // cout  << "Center: " << center << endl;
  // cout  << "Sorted Lattice IDs: "; print1DVector(sorted_lattice_ids);

  return sorted_lattice_ids;
}

struct Vector3dEqual
{
  bool operator()(const Vector3d &a, const Vector3d &b) const
  {
    return (a - b).norm() < K_EPSILON;
  }
};

vector<vector<size_t>> GetEquivalentLatticeIdUnder3BarSymmetry(
    const Config &config,
    size_t maxBondOrder,
    const pair<size_t, size_t> &latticeIdPair)
{
  // Get sorted lattice IDs
  auto ssVector = config.GetSortedLatticeVectorStateOfPair(latticeIdPair, maxBondOrder);

  // Compute pair direction and rotation axis
  Vector3d centralPos = config.GetRelativePositionOfLattice(latticeIdPair.first);
  Vector3d nnPos = config.GetRelativePositionOfLattice(latticeIdPair.second);
  Vector3d pairDirection = centralPos - nnPos;
  Vector3d rotationAxis = GetClosest111Direction(pairDirection);

  // Compute center
  Vector3d center = GetLatticePairCenter(config, latticeIdPair).transpose();

  // Store positions and indices
  unordered_map<Vector3d, size_t, Vector3dHash> positionToIndex;
  vector<Vector3d> relativePositionVector;
  relativePositionVector.reserve(ssVector.size());

  for (size_t i = 0; i < ssVector.size(); ++i)
  {
    auto pos = config.GetRelativePositionOfLattice(ssVector[i]);
    relativePositionVector.push_back(pos);
    positionToIndex[pos] = i;
  }

  // Track processed sites
  vector<bool> processed(ssVector.size(), false);
  vector<vector<size_t>> equivalentEncodingVector;
  vector<vector<size_t>> equivalentLatticeIdVector;

  for (size_t i = 0; i < ssVector.size(); ++i)
  {
    if (processed[i])
    {
      continue;
    }

    vector<size_t> equivalentIndices = {i};
    vector<size_t> equivalentLatticeIds = {ssVector[i]};
    processed[i] = true;

    // Compute equivalent positions with PBC
    auto equivPositions = GetEquivalentPositionUnder3BarSymmetry(
        relativePositionVector[i], rotationAxis, center);

    for (const auto &equivPos : equivPositions)
    {
      size_t index = findClosestMatch(positionToIndex, equivPos);
      if (index != numeric_limits<size_t>::max() && !processed[index])
      {
        equivalentIndices.push_back(index);
        equivalentLatticeIds.push_back(ssVector[index]);

        processed[index] = true;
      }
    }

    // Sort indices within each equivalent set
    // sort(equivalentIndices.begin(), equivalentIndices.end());
    equivalentEncodingVector.push_back(move(equivalentIndices));
    equivalentLatticeIdVector.push_back(move(equivalentLatticeIds));
  }

  // Sort equivalent sets by the first index
  // sort(equivalentEncodingVector.begin(), equivalentEncodingVector.end(),
  //      [](const auto &a, const auto &b)
  //      { return a[0] < b[0]; });

  print2DVector(equivalentEncodingVector);
  print2DVector(equivalentLatticeIdVector);

  return equivalentLatticeIdVector;
}

////////////////////////////////////

////// Custom Function ////////////

// void wrap(Vector3d &relativePosition)
// {
//   for (const int kDim : std::vector<int>{0, 1, 2})
//   {
//     while (relativePosition[kDim] >= 1)
//     {
//       relativePosition[kDim] -= 1;
//     }
//     while (relativePosition[kDim] < 0)
//     {
//       relativePosition[kDim] += 1;
//     }
//   }
// }

template <typename Derived>
void wrap(Eigen::MatrixBase<Derived> &vec)
{
  vec = vec.unaryExpr([](double x)
                      { return x - std::floor(x); });
}

Vector3d GetEquivalentPositionUnderInversion(const Vector3d &position,
                                             const Vector3d &inversionCenter)
{
  Vector3d invertedPosition = 2 * inversionCenter - position;
  return invertedPosition;
}

bool areApproximatelyEqual(const Vector3d &v1,
                           const Vector3d &v2)
{
  return (v1 - v2).cwiseAbs().maxCoeff() < constants::kEpsilon;
}

size_t GetLatticeIdByPosition(
    const Vector3d &position,
    const unordered_map<Vector3d, size_t, Vector3dHash> &positionLatticeIdMap)
{
  for (const auto &entry : positionLatticeIdMap)
  {
    if (areApproximatelyEqual(position, entry.first))
    {
      return entry.second;
    }
  }

  // Need to work on this
  cerr << "No equivalent site found!" << endl;
  // exit(3);
  // throw runtime_error("Equivalent position not found in lattice.");
}

set<pair<size_t, size_t>> GetEquivalentSitesByInversion(
    const Vector3d &inversionCenter,
    const unordered_map<Vector3d, size_t, Vector3dHash> &positionToLatticeIdMap)
{
  set<pair<size_t, size_t>> equivalentSitePairs;

  for (const auto &entry : positionToLatticeIdMap)
  {
    const Vector3d &originalPosition = entry.first;
    size_t originalId = entry.second;

    Vector3d invertedPosition = GetEquivalentPositionUnderInversion(originalPosition,
                                                                    inversionCenter);

    size_t invertedId = GetLatticeIdByPosition(invertedPosition,
                                               positionToLatticeIdMap);

    size_t id1 = originalId;
    size_t id2 = invertedId;
    if (id1 > id2)
      std::swap(id1, id2);
    equivalentSitePairs.emplace(id1, id2);
  }

  return equivalentSitePairs;
}

Matrix3d GetRotationMatrix(const Vector3d &rotationAxis,
                           const double rotationAngle)
{
  // Using rodrigues formula
  Vector3d normAxis = rotationAxis.normalized();

  double ux = normAxis.x(), uy = normAxis.y(), uz = normAxis.z();

  auto cosTheta = cos(toRadians(rotationAngle));
  auto sinTheta = sin(toRadians(rotationAngle));

  Matrix3d rotationMatrix;

  rotationMatrix(0, 0) = cosTheta + ux * ux * (1 - cosTheta);
  rotationMatrix(0, 1) = ux * uy * (1 - cosTheta) - uz * sinTheta;
  rotationMatrix(0, 2) = ux * uz * (1 - cosTheta) + uy * sinTheta;

  rotationMatrix(1, 0) = uy * ux * (1 - cosTheta) + uz * sinTheta;
  rotationMatrix(1, 1) = cosTheta + uy * uy * (1 - cosTheta);
  rotationMatrix(1, 2) = uy * uz * (1 - cosTheta) - ux * sinTheta;

  rotationMatrix(2, 0) = uz * ux * (1 - cosTheta) - uy * sinTheta;
  rotationMatrix(2, 1) = uz * uy * (1 - cosTheta) + ux * sinTheta;
  rotationMatrix(2, 2) = cosTheta + uz * uz * (1 - cosTheta);

  return rotationMatrix;
}

set<pair<size_t, size_t>> GetEquivalentSitesUnderRotation(
    const Vector3d &rotationAxis,
    const unordered_map<Vector3d, size_t, Vector3dHash> &positionToLatticeIdMap)
{
  // Rotation angle
  const double rotationAngle = 120; // degree

  // Rotation matrix
  Matrix3d rotationMatrix = GetRotationMatrix(rotationAxis,
                                              rotationAngle);

  cout << "Rotation Matrix" << endl;
  cout << rotationMatrix << endl;

  set<pair<size_t, size_t>> equivalentSitePairs;

  for (const auto &entry : positionToLatticeIdMap)
  {
    const Vector3d &originalPosition = entry.first;
    size_t originalId = entry.second;

    Vector3d rotatedPosition = rotationMatrix * originalPosition;
    wrap(rotatedPosition);

    size_t rotatedId = GetLatticeIdByPosition(rotatedPosition,
                                              positionToLatticeIdMap);

    size_t id1 = originalId;
    size_t id2 = rotatedId;
    if (id1 > id2)
      std::swap(id1, id2);
    equivalentSitePairs.emplace(id1, id2);
  }

  return equivalentSitePairs;
}

#include "UnionFind.h"

Vector3d GetLatticePairCenterCorrect(
    const Config &config,
    const std::pair<size_t, size_t> &latticeIdPair)
{

  // Get relative positions of the lattice IDs
  Eigen::Vector3d position1 = config.GetRelativePositionOfLattice(latticeIdPair.first);
  Eigen::Vector3d position2 = config.GetRelativePositionOfLattice(latticeIdPair.second);

  Eigen::Vector3d centerPosition;
  for (int kDim = 0; kDim < 3; ++kDim)
  {
    double distance = position1[kDim] - position2[kDim];
    int period = static_cast<int>(distance / 0.5);

    // Adjust the positions once based on the period
    if (period != 0)
    {
      position1[kDim] -= period;
    }

    centerPosition[kDim] = 0.5 * (position1[kDim] + position2[kDim]);
  }

  return centerPosition;
}

#include <vector>
#include <array>
#include <algorithm>
#include <Eigen/Dense>
// Helper function to compute signed tetrahedron volume
double signed_tetra_volume(const Eigen::Vector3d &a,
                           const Eigen::Vector3d &b,
                           const Eigen::Vector3d &c,
                           const Eigen::Vector3d &d)
{
  return (a - d).dot((b - d).cross(c - d)) / 6.0;
}
double convex_hull_volume_6pts(const std::array<Eigen::Vector3d, 6> &points)
{
  // Step 1: Find the point with minimum x-coordinate as our base
  auto base_it = std::min_element(points.begin(), points.end(),
                                  [](const Eigen::Vector3d &a, const Eigen::Vector3d &b)
                                  {
                                    return a.x() < b.x();
                                  });
  const Eigen::Vector3d base = *base_it;
  size_t base_index = std::distance(points.begin(), base_it);

  // Step 2: Compute volumes of all possible tetrahedrons containing the base point
  std::array<double, 10> volumes;
  size_t idx = 0;

  for (size_t i = 0; i < 6; ++i)
  {
    if (i == base_index)
      continue;
    for (size_t j = i + 1; j < 6; ++j)
    {
      if (j == base_index)
        continue;
      for (size_t k = j + 1; k < 6; ++k)
      {
        if (k == base_index)
          continue;
        double vol = signed_tetra_volume(base, points[i], points[j], points[k]);
        if (vol > 0) // Only consider positive volume for outer tetrahedrons
          volumes[idx++] = vol;
      }
    }
  }

  // Step 3: Sum the positive volumes
  double total_volume = 0.0;
  for (size_t i = 0; i < idx; ++i)
  {
    total_volume += volumes[i];
  }

  return total_volume;
}

#include <Eigen/Eigenvalues>

// Helper function to compute moment of inertia tensor eigenvalues
double get_inertia_eigenvalue_sum(const std::array<Eigen::Vector3d, 6> &positions)
{
  // Step 1: Compute center of mass
  Eigen::Vector3d com(0.0, 0.0, 0.0);
  for (const auto &p : positions)
  {
    com += p;
  }
  com /= 6.0;

  cout << "COM : " << com.transpose() << endl;

  // Step 2: Compute moment of inertia tensor
  Eigen::Matrix3d I = Eigen::Matrix3d::Zero();
  for (const auto &p : positions)
  {
    Eigen::Vector3d r = p - com;
    I(0, 0) += r(1) * r(1) + r(2) * r(2);
    I(1, 1) += r(0) * r(0) + r(2) * r(2);
    I(2, 2) += r(0) * r(0) + r(1) * r(1);
    I(0, 1) -= r(0) * r(1);
    I(0, 2) -= r(0) * r(2);
    I(1, 2) -= r(1) * r(2);
  }
  I(1, 0) = I(0, 1);
  I(2, 0) = I(0, 2);
  I(2, 1) = I(1, 2);

  // Step 3: Compute eigenvalues
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(I);
  Eigen::Vector3d eigenvalues = solver.eigenvalues();

  // Step 4: Return sum of eigenvalues (or sort them for consistency)
  return eigenvalues.sum(); // Alternatively, use eigenvalues(0) or a sorted combination
}

void writeToDataFile(const pair<size_t, size_t> &latticeIdJumpPair,
                     const unordered_map<Vector3d, size_t, Vector3dHash> &positionLatticeIdMap,
                     vector<vector<size_t>> &equivalentLatticeIdVector)
// const vector<pair<double, vector<size_t>>> &sortable_classes)
{
  ofstream outFile;
  string filename = to_string(latticeIdJumpPair.first) + "_" +
                    to_string(latticeIdJumpPair.second) + "_reducedPositionLCE.dump";

  outFile.open(filename);

  outFile << "ITEM: TIMESTEP\n0\n";
  outFile << "ITEM: NUMBER OF ATOMS\n"
          << positionLatticeIdMap.size() << "\n";
  outFile << "ITEM: BOX BOUNDS pp pp pp\n0 1\n0 1\n0 1\n";
  outFile << "ITEM: ATOMS id type x y z group\n";

  // Create a lookup table for quick ID-to-class_value mapping
  unordered_map<size_t, int> idToClassValue;
  int i = 0;
  for (const auto &eqSites : equivalentLatticeIdVector)
  {
    int group = i;
    for (size_t site : eqSites)
    {
      idToClassValue[site] = group;
    }
    i += 1;
  }

  for (const auto &entry : positionLatticeIdMap)
  {
    size_t id = entry.second;
    outFile << id << " 1 ";
    for (const auto &coord : entry.first)
    {
      outFile << coord << " ";
    }

    // Add the class value or -1 if not found
    if (idToClassValue.count(id))
    {
      outFile << idToClassValue[id];
    }
    else
    {
      outFile << -1.0; // or some default/error value
    }

    outFile << "\n";
  }

  outFile.close();
}
/*
#include <qhull/qhull.h>
#include <qhull/qhull_r.h>


// Function to compute convex hull volume using Qhull
double compute_convex_hull_volume(const std::array<Eigen::Vector3d, 6>& positions) {
    // Convert Eigen vectors to Qhull points
    std::vector<double> points;
    for (const auto& pos : positions) {
        points.push_back(pos.x());
        points.push_back(pos.y());
        points.push_back(pos.z());
    }

    // Set up Qhull for 3D points
    Qhull qhull;
    QhullError err;

    // Run Qhull to compute the convex hull
    qhull.runQhull("", false, false, false, points, err);

    if (err.isError()) {
        std::cerr << "Qhull error: " << err.getMessage() << std::endl;
        return 0.0;
    }

    // Get the volume of the convex hull
    double volume = qhull.getVolume();

    return volume;
}
    */

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <utility>
#include <unordered_map>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <utility>
#include <unordered_map>

vector<vector<size_t>> GetEquivalentSitesUnder3BarSymmetry(const Config &config,
                                              const pair<size_t, size_t> latticeIdJumpPair,
                                              const size_t maxBondOrder)
{
  auto ssVector = config.GetSortedLatticeVectorStateOfPair(latticeIdJumpPair,
                                                           maxBondOrder);

  // Position to lattice Id map
  unordered_map<Vector3d, size_t, Vector3dHash> positionLatticeIdMap;

  unordered_map<size_t, Vector3d> idPositionMap;

  vector<Vector3d> modifiedPositionVector;

  // Move all the position to center
  // Vector3d moveToCenter(0.5, 0.5, 0.5);
  Vector3d position1 = config.GetRelativePositionOfLattice(latticeIdJumpPair.first);
  Vector3d position2 = config.GetRelativePositionOfLattice(latticeIdJumpPair.second);

  Vector3d moveToCenter = Vector3d(0.5, 0.5, 0.5) - GetLatticePairCenterCorrect(config, latticeIdJumpPair);

  // wrap(inversionCenter);

  Vector3d modifiedPosition1 = position1 + moveToCenter;
  Vector3d modifiedPosition2 = position2 + moveToCenter;

  wrap(modifiedPosition1);
  wrap(modifiedPosition2);

  Vector3d inversionCenter = (modifiedPosition1 + modifiedPosition2) / 2.0; // {0.5, 0.5, 0.5}

  // Rotation axis
  Vector3d rotationAxis = modifiedPosition1 - modifiedPosition2;
  // wrap(rotationAxis);

  cout << "Rotation Axis : " << rotationAxis.transpose() << endl;

  for (const auto id : ssVector)
  {
    Vector3d modifiedPosition = config.GetRelativePositionOfLattice(id);
    modifiedPosition += moveToCenter; // Shift to center (0.5, 0.5, 0.5)
    wrap(modifiedPosition);

    modifiedPositionVector.emplace_back(modifiedPosition);

    positionLatticeIdMap[modifiedPosition] = id;

    // test
    idPositionMap[id] = modifiedPosition; // Keep this updated too

    // cout << id << " : " << modifiedPosition.transpose() << endl;
  }

  // Equivalent sites under inversion
  auto equivalentSitesInversion = GetEquivalentSitesByInversion(inversionCenter,
                                                                positionLatticeIdMap);

  vector<vector<size_t>> eqSites3BarVector;

  // cout << "Equivalent Sites Under Inversion" << endl;
  for (auto eqSite : equivalentSitesInversion)
  {
    // cout << eqSite.first << ", " << eqSite.second << endl;
    eqSites3BarVector.emplace_back(vector<size_t>{eqSite.first, eqSite.second});
  }

  // Equivalent sites under Rotation

  auto equivalentSitesRotation = GetEquivalentSitesUnderRotation(rotationAxis,
                                                                 positionLatticeIdMap);

  // vector<vector<size_t>> equivalentSitesRotationVector;

  // cout << "Equivalent Sites Under Rotation" << endl;
  for (auto eqSite : equivalentSitesRotation)
  {
    // cout << eqSite.first << ", " << eqSite.second << endl;
    eqSites3BarVector.emplace_back(vector<size_t>{eqSite.first, eqSite.second});
  }

  // Combined equivalent site under inversion and rotation

  // Combine the vectors

  // print2DVector(equivalentSitesRotationVector);
  // print2DVector(equivalentSitesInversionVector);

  // Find maximum index
  UnionFind<size_t> uf;

  // Process combined relations
  for (const auto &eqSites : eqSites3BarVector)
  {
    if (eqSites.size() == 2)
    {
      uf.union_sets(eqSites[0], eqSites[1]);
    }
    else if (eqSites.size() == 1)
    {
      uf.find(eqSites[0]); // Ensure it's in the parent map
    }
  }

  // First, build map from root to members
  unordered_map<size_t, vector<size_t>> eqSites3BarMap;

  auto equivalenceSiteMap = uf.getEquivalenceMap();

  for (const auto &entry : equivalenceSiteMap)
  {
    auto root = uf.find(entry.first);
    eqSites3BarMap[root].push_back(entry.first);
  }

  // Convert to vector of vector
  vector<vector<size_t>> equivalence_classes;
  for (auto &kv : eqSites3BarMap)
  {
    sort(kv.second.begin(), kv.second.end());
    equivalence_classes.push_back(kv.second);
  }

  // Using Volume for sorting
  cout << "Equivalence Lattice Ids" << endl;

//  std::vector<std::pair<double, std::vector<size_t>>> sortable_classes;
//
//  for (auto &eq_class : equivalence_classes)
//  {
//    if (eq_class.size() == 6)
//    {
//      std::array<Eigen::Vector3d, 6> positions;
//      for (size_t i = 0; i < 6; ++i)
//      {
//        positions[i] = idPositionMap.at(eq_class[i]);
//      }
//
//      double vol = convex_hull_volume_6pts(positions);
//      sortable_classes.emplace_back(vol, eq_class);
//    }
//    else if (eq_class.size() == 2)
//    {
//      sortable_classes.emplace_back(0, eq_class);
//    }
//  }
//
//  // Step 2: Sort by volume
//  std::sort(sortable_classes.begin(), sortable_classes.end(),
//            [](const auto &a, const auto &b)
//            {
//              return a.first < b.first; // ascending order
//            });
//
//  // Step 3: Extract sorted equivalence classes


  // cout << "Equivalent Lattice Ids" << endl;
  // print2DVector(sorted_equivalence_classes);
  // print2DVector(equivalentSites);

  // Modified sorting section in GetEquivalentSites3Bar
  std::vector<std::pair<double, std::vector<size_t>>> sortable_classes;
  for (auto &eq_class : equivalence_classes)
  {
    if (eq_class.size() == 6)
    {
      std::array<Eigen::Vector3d, 6> positions;
      for (size_t i = 0; i < 6; ++i)
      {
        positions[i] = idPositionMap.at(eq_class[i]);
      }
      double inertia_sum = get_inertia_eigenvalue_sum(positions);
      sortable_classes.emplace_back(inertia_sum, eq_class);
    }
    else if (eq_class.size() == 2)
    {
      sortable_classes.emplace_back(0, eq_class); // Keep pairs at the start
    }
  }

  // Sort by inertia eigenvalue sum
  std::sort(sortable_classes.begin(), sortable_classes.end(),
            [](const auto &a, const auto &b)
            {
              return a.first < b.first; // Ascending order
            });

    std::vector<std::vector<size_t>> sorted_equivalence_classes;
    for (const auto &entry : sortable_classes)
    {
      sorted_equivalence_classes.push_back(entry.second);
  
      cout << " { ";
      for (auto id : entry.second)
      {
        cout << id << ", ";
      }
      cout << " } : " << entry.first << endl;
    }

  // writeToDataFile(latticeIdJumpPair,
  //                 positionLatticeIdMap,
  //                 sorted_equivalence_classes);

  return sorted_equivalence_classes;
}

// #include <libqhull/geom.h> // Removed as it is not found