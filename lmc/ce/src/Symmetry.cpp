#include "Symmetry.h"

// Convert degrees to radians
inline double toRadians(double degrees)
{
  return degrees * M_PI / 180.0;
}

// Wraps the position to [0,1)
template <typename Derived>
void wrap(Eigen::MatrixBase<Derived> &vec)
{
  vec = vec.unaryExpr([](double x)
                      { return x - floor(x); });
}

// Get equivalent sites under inversion
Vector3d GetEquivalentPositionUnderInversion(const Vector3d &position,
                                             const Vector3d &inversionCenter)
{
  Vector3d invertedPosition = 2 * inversionCenter - position;
  return invertedPosition;
}

// Compares two vectors and return true if the vectors are equal
bool areApproximatelyEqual(const Vector3d &v1,
                           const Vector3d &v2)
{
  return (v1 - v2).cwiseAbs().maxCoeff() < constants::kEpsilon;
}

// Returns the lattice ID given the position
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
  exit(3);
  // throw runtime_error("Equivalent position not found in lattice.");
}

// Maps the sites if they are equivalent under inversion
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
      swap(id1, id2);
    equivalentSitePairs.emplace(id1, id2);
  }

  return equivalentSitePairs;
}

// Returns the rotation matrix using rodrigues formula
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

// Maps the sites if they are equivalent under rotation
set<pair<size_t, size_t>> GetEquivalentSitesUnderRotation(
    const Vector3d &rotationAxis,
    const double rotationAngle,
    const unordered_map<Vector3d, size_t, Vector3dHash> &positionToLatticeIdMap)
{
  // Rotation matrix
  Matrix3d rotationMatrix = GetRotationMatrix(rotationAxis,
                                              rotationAngle);

  // cout << "Rotation Matrix" << endl;
  // cout << rotationMatrix << endl;

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
      swap(id1, id2);
    equivalentSitePairs.emplace(id1, id2);
  }

  return equivalentSitePairs;
}

// Lattice pair center
Vector3d GetLatticePairCenterCorrect(
    const Config &config,
    const pair<size_t, size_t> &latticeIdPair)
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

// Function to compute moment of inertia tensor eigenvalues
// Used for sorting the equivalent classes which will be rotation invariant
double GetInteriaEigenValueSum(const array<Eigen::Vector3d, 6> &positions)
{
  // Step 1: Compute center of mass
  Eigen::Vector3d com(0.0, 0.0, 0.0);
  for (const auto &p : positions)
  {
    com += p;
  }
  com /= 6.0;

  // cout << "COM : " << com.transpose() << endl;

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
  // Alternatively, use eigenvalues(0) or a sorted combination
  return eigenvalues.sum();
}

// Returns equivalent sites under 3 Bar symmetry
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

  Vector3d modifiedPosition1 = position1 + moveToCenter;
  Vector3d modifiedPosition2 = position2 + moveToCenter;

  wrap(modifiedPosition1);
  wrap(modifiedPosition2);

  Vector3d inversionCenter = (modifiedPosition1 + modifiedPosition2) / 2.0; // {0.5, 0.5, 0.5}

  // Rotation axis
  Vector3d rotationAxis = modifiedPosition1 - modifiedPosition2;

  // cout << "Rotation Axis : " << rotationAxis.transpose() << endl;

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

  // Rotation angle
  const double rotationAngle = 120; // degree
  auto equivalentSitesRotation = GetEquivalentSitesUnderRotation(rotationAxis,
                                                                 rotationAngle,
                                                                 positionLatticeIdMap);

  for (auto eqSite : equivalentSitesRotation)
  {
    eqSites3BarVector.emplace_back(vector<size_t>{eqSite.first, eqSite.second});
  }

  // Combined equivalent site under inversion and rotation

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
      uf.find(eqSites[0]);
    }
  }

  // Merging the equivalent sites from 3 fold rotation and inversion
  unordered_map<size_t, vector<size_t>> eqSites3BarMap;

  auto equivalenceSiteMap = uf.getEquivalenceMap();

  for (const auto &entry : equivalenceSiteMap)
  {
    auto root = uf.find(entry.first);
    eqSites3BarMap[root].push_back(entry.first);
  }

  // Convert to vector of vector
  vector<vector<size_t>> equivalentClasses3Bar;
  for (auto &kv : eqSites3BarMap)
  {
    sort(kv.second.begin(), kv.second.end());
    equivalentClasses3Bar.push_back(kv.second);
  }

  // cout << "Equivalence Lattice Ids" << endl;

  // Sorting the equivalent site using interia eigen value sum
  vector<pair<double, vector<size_t>>> sortableEquivalentClasses;

  for (auto &eqClass : equivalentClasses3Bar)
  {
    if (eqClass.size() == 6)
    {
      array<Eigen::Vector3d, 6> positions;
      for (size_t i = 0; i < 6; ++i)
      {
        positions[i] = idPositionMap.at(eqClass[i]);
      }
      // Inertia sum
      double inertiaSum = GetInteriaEigenValueSum(positions);
      sortableEquivalentClasses.emplace_back(inertiaSum, eqClass);
    }
    else if (eqClass.size() == 2)
    {
      sortableEquivalentClasses.emplace_back(0, eqClass);
    }
  }

  // Sort by inertia eigenvalue sum
  sort(sortableEquivalentClasses.begin(), sortableEquivalentClasses.end(),
       [](const auto &a, const auto &b)
       {
         // Ascending order
         return a.first < b.first;
       });

  vector<vector<size_t>> sortedEquivalentClass3Bar;
  for (const auto &entry : sortableEquivalentClasses)
  {
    sortedEquivalentClass3Bar.push_back(entry.second);

    // cout << " { ";
    // for (auto id : entry.second)`
    // {
    //   cout << id << ", ";
    // }
    // cout << " } : " << entry.first << endl;
  }

  return sortedEquivalentClass3Bar;
}

// Compute the rotation matrix to align `currentAxis` with `targetAxis`
Matrix3d GetRotationMatrixToAlignVectors(
    const Vector3d &currentAxis,
    const Vector3d &targetAxis)
{
  Vector3d v1 = currentAxis.normalized();
  Vector3d v2 = targetAxis.normalized();
  Vector3d v = v1.cross(v2);
  double v_norm = v.norm();

  // Handle special case: vectors are parallel or anti-parallel
  if (v_norm < constants::kEpsilon)
  {
    double dot = v1.dot(v2);
    if (dot > 0.0)
    {
      return Matrix3d::Identity(); // No rotation needed
    }
    else
    {
      // 180-degree rotation around a perpendicular axis
      Vector3d perp;
      if (std::abs(v1.x()) > std::abs(v1.y()))
      {
        perp = Vector3d(-v1.y(), v1.x(), 0);
      }
      else
      {
        perp = Vector3d(0, -v1.z(), v1.y());
      }
      perp.normalize();
      return Matrix3d::Identity() - 2.0 * perp * perp.transpose();
    }
  }

  // Normalize rotation axis
  v.normalize();

  double cos_theta = v1.dot(v2);
  cos_theta = std::clamp(cos_theta, -1.0, 1.0); // Clamp to avoid numerical errors
  double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);

  // Skew-symmetric matrix for cross product
  Matrix3d K;
  K << 0, -v.z(), v.y(),
      v.z(), 0, -v.x(),
      -v.y(), v.x(), 0;

  // Rodrigues' rotation formula
  return Matrix3d::Identity() + sin_theta * K + (1 - cos_theta) * K * K;
}

vector<size_t> GetSSVector3FSymmetry(
    const Config &config,
    const pair<size_t, size_t> &latticeIdJumpPair,
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

  Vector3d modifiedPosition1 = position1 + moveToCenter;
  Vector3d modifiedPosition2 = position2 + moveToCenter;

  wrap(modifiedPosition1);
  wrap(modifiedPosition2);

  Vector3d inversionCenter = (modifiedPosition1 + modifiedPosition2) / 2.0; // {0.5, 0.5, 0.5}

  // Rotation axis
  Vector3d rotationAxis = modifiedPosition1 - modifiedPosition2;

  // cout << "Rotation Axis : " << rotationAxis.transpose() << endl;

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

  Vector3d targetAxis(0, 0, 1);

  Matrix3d rotationMatrix = GetRotationMatrixToAlignVectors(rotationAxis,
                                                            targetAxis);

  vector<tuple<size_t, Vector3d, Vector3d, Vector3d>> transformedEntries;

  for (const auto &entry : positionLatticeIdMap)
  {
    Vector3d original = entry.first;
    Vector3d reduced = original - Vector3d(0.5, 0.5, 0.5);
    Vector3d transformed = rotationMatrix * reduced;

    transformedEntries.emplace_back(entry.second, original, reduced, transformed);
  }

  // Sort based on transformed position: first z, then y, then x
  sort(transformedEntries.begin(), transformedEntries.end(),
       [](const auto &a, const auto &b)
       {
         const Vector3d &ta = get<3>(a);
         const Vector3d &tb = get<3>(b);
         if (ta.z() != tb.z())
           return ta.z() < tb.z();
         if (ta.y() != tb.y())
           return ta.y() < tb.y();
         return ta.x() < tb.x();
       });

  // Print
  // int eId = 0;
  // cout << "Id\teId\tOriginal\tReduced\tTransformed" << endl;
  // for (const auto &entry : transformedEntries)
  // {
  //   cout << get<0>(entry) << "\t"
  //        << eId << "\t"
  //        << get<1>(entry).transpose() << "\t"
  //        << get<2>(entry).transpose() << "\t"
  //        << get<3>(entry).transpose() << endl;
  //   eId++;
  // }
  vector<size_t> ssVector3Bar;
  for (const auto &entry : transformedEntries)
  {
    ssVector3Bar.emplace_back(get<0>(entry));
  }

  return ssVector3Bar;
}

vector<vector<size_t>> GetEquivalentSiteEncoding3BarSymmetry(
    const Config &config,
    size_t maxBondOrder)
{
  // Get the equivalent lattice Ids
  pair<size_t, size_t> latticeIdJumpPair = {0, config.GetNeighborLatticeIdVectorOfLattice(0, 1)[0]};

  auto equivalentLatticeIds3Bar = GetEquivalentSitesUnder3BarSymmetry(
      config,
      latticeIdJumpPair,
      maxBondOrder);

  auto ssVector = GetSSVector3FSymmetry(
      config,
      latticeIdJumpPair,
      maxBondOrder);

  unordered_map<size_t, size_t> latticeToIndexMap;

  size_t idx = 0;
  for (const size_t &latticeId : ssVector)
  {
    latticeToIndexMap[latticeId] = idx;
    idx++;
  }

  // Equivalent sites encoding
  vector<vector<size_t>> equivalentSitesEncoding3Bar;
  equivalentSitesEncoding3Bar.reserve(equivalentLatticeIds3Bar.size());

  for (const auto &eqGroup : equivalentLatticeIds3Bar)
  {
    vector<size_t> localEncoding;
    localEncoding.reserve(eqGroup.size());

    for (auto const &latticeId : eqGroup)
    {
      localEncoding.emplace_back(latticeToIndexMap.at(latticeId));
    }
    sort(localEncoding.begin(), localEncoding.end());
    equivalentSitesEncoding3Bar.emplace_back(localEncoding);
  }

  return equivalentSitesEncoding3Bar;
}