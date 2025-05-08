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
 * @brief Computes the cartesian center position of a lattice pair in a periodic lattice system.
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

  return (center_position.transpose() * config.GetBasis()).transpose();
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
