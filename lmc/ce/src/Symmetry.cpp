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
    if (dot > max_dot + constants::kEpsilon)
    {
      max_dot = dot;
      best_direction = dir;
    }
    else if (fabs(dot - max_dot) < constants::kEpsilon)
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
static Vector3d wrapToUnitCell(const Vector3d &pos)
{
  Vector3d wrapped = pos;
  for (int i = 0; i < 3; ++i)
  {
    wrapped[i] = fmod(wrapped[i], 1.0);
    if (wrapped[i] < 0.0)
    {
      wrapped[i] += 1.0;
    }
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
 * within a predefined tolerance (constants::kEpsilon), the corresponding index is returned.
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
    if (periodicDistance(pos, target) < constants::kEpsilon)
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
      if (a[i] + constants::kEpsilon < b[i]) return true;
      if (a[i] - constants::kEpsilon > b[i]) return false;
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
    if (diff < -constants::kEpsilon)
      return true;
    if (diff > constants::kEpsilon)
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
  Config::RotateLatticeVector(lattice_id_hashmap, config.GetLatticePairRotationMatrix(lattice_id_jump_pair));

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
