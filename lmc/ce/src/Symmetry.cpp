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

  //cout << "Equivalence Lattice Ids" << endl;
  

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
