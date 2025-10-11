#ifndef LMC_CFG_INCLUDE_CUBE_H_
#define LMC_CFG_INCLUDE_CUBE_H_

#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

class Cube
{
public:
  Cube(
      const size_t &sizeOfCube);

  const vector<size_t> &GetNeighbors(size_t siteIndex) const;

  Vector3i GetRelativePosition(size_t siteIndex) const;

  Matrix3i GetBasis() const;

  Matrix3Xi GetRelativePositionMatrix() const;

  size_t GetSizeOfCube() const;

  size_t GetNumOfSites() const;

  size_t GetCentralSiteId() const;

private:
  void UpdateNeighbors();

  const size_t sizeOfCube_{};
  Matrix3i basis_{};                // factor
  Matrix3Xi relativeIndexMatrix_{}; // will store the (i, j, k)
  vector<vector<size_t>> neighbors_;
};

#endif // LMC_CFG_INCLUDE_CUBE_H_
