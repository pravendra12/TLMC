#include "Cube.h"
#include "Constants.hpp"

Cube::Cube(const size_t &sizeOfCube) : sizeOfCube_(sizeOfCube)
{
  basis_ << sizeOfCube_, 0, 0,
      0, sizeOfCube_, 0,
      0, 0, sizeOfCube_;

  size_t totalCubes = sizeOfCube_ * sizeOfCube_ * sizeOfCube_;
  relativeIndexMatrix_.resize(3, static_cast<int>(totalCubes));

  size_t colIdx = 0;
  for (int ii = 0; ii < sizeOfCube_; ++ii)
  {
    for (int jj = 0; jj < sizeOfCube_; ++jj)
    {
      for (int kk = 0; kk < sizeOfCube_; ++kk)
      {
        relativeIndexMatrix_.col(colIdx) = Vector3i(ii, jj, kk);
        colIdx++;
      }
    }
  }

  UpdateNeighbors();
}
/*
void Cube::UpdateNeighbors()
{
  size_t totalCubes = sizeOfCube_ * sizeOfCube_ * sizeOfCube_;
  neighbors_.clear();
  neighbors_.resize(totalCubes);

  for (size_t idx = 0; idx < totalCubes; ++idx)
  {
    Vector3i pos = relativeIndexMatrix_.col(idx);

    // Loop over neighbor offsets in a 3x3x3 cube (excluding self)
    for (int dx = -1; dx <= 1; ++dx)
    {
      for (int dy = -1; dy <= 1; ++dy)
      {
        for (int dz = -1; dz <= 1; ++dz)
        {
          if (dx == 0 && dy == 0 && dz == 0)
            continue; // skip self

          // Apply periodic boundary conditions
          int nx = (pos(0) + dx + sizeOfCube_) % sizeOfCube_;
          int ny = (pos(1) + dy + sizeOfCube_) % sizeOfCube_;
          int nz = (pos(2) + dz + sizeOfCube_) % sizeOfCube_;

          size_t neighborIdx = nx * sizeOfCube_ * sizeOfCube_ +
                               ny * sizeOfCube_ +
                               nz;
          neighbors_[idx].push_back(neighborIdx);
        }
      }
    }
  }
}
  */

inline bool PositionCompareState(const pair<size_t, Eigen::RowVector3d> &lhs,
                                 const pair<size_t, Eigen::RowVector3d> &rhs)
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

void Cube::UpdateNeighbors()
{
  size_t totalCubes = sizeOfCube_ * sizeOfCube_ * sizeOfCube_;
  neighbors_.clear();
  neighbors_.resize(totalCubes);

  for (size_t idx = 0; idx < totalCubes; ++idx)
  {
    Vector3i pos = relativeIndexMatrix_.col(idx);

    // --- Step 1: Store neighbor IDs ---
    vector<size_t> neighborIds;
    for (int dx = -1; dx <= 1; ++dx)
    {
      for (int dy = -1; dy <= 1; ++dy)
      {
        for (int dz = -1; dz <= 1; ++dz)
        {
          if (dx == 0 && dy == 0 && dz == 0)
            continue;

          int nx = (pos(0) + dx + sizeOfCube_) % sizeOfCube_;
          int ny = (pos(1) + dy + sizeOfCube_) % sizeOfCube_;
          int nz = (pos(2) + dz + sizeOfCube_) % sizeOfCube_;

          size_t neighborIdx = nx * sizeOfCube_ * sizeOfCube_ +
                               ny * sizeOfCube_ +
                               nz;

          neighborIds.push_back(neighborIdx);
        }
      }
    }

    // --- Step 2: Center neighbor positions ---
    vector<pair<size_t, Eigen::RowVector3d>> centeredNeighbors;
    centeredNeighbors.reserve(neighborIds.size());

    Eigen::RowVector3d centerPos(
        double(pos(0)) / sizeOfCube_,
        double(pos(1)) / sizeOfCube_,
        double(pos(2)) / sizeOfCube_);

    for (size_t nId : neighborIds)
    {
      Vector3i nPos = relativeIndexMatrix_.col(nId);

      Eigen::RowVector3d relPos(
          double(nPos(0)) / sizeOfCube_,
          double(nPos(1)) / sizeOfCube_,
          double(nPos(2)) / sizeOfCube_);

      // Center around current site
      relPos -= centerPos;

      // Minimal image symmetric [-0.5,0.5)
      relPos = relPos.unaryExpr([](double x)
                                { return x - std::floor(x + 0.5); });

      // Shift to [0,1)
      relPos += Eigen::RowVector3d(0.5, 0.5, 0.5);

      centeredNeighbors.emplace_back(nId, relPos);
    }

    // --- Step 3: Sort neighbors canonically ---
    std::sort(centeredNeighbors.begin(), centeredNeighbors.end(),PositionCompareState);

    // Store sorted neighbor IDs
    neighbors_[idx].reserve(centeredNeighbors.size());
    for (const auto &p : centeredNeighbors)
      neighbors_[idx].push_back(p.first);
  }
}

const vector<size_t>& Cube::GetNeighbors(size_t siteIndex) const
{
  if (siteIndex >= neighbors_.size())
    throw std::out_of_range("Site index out of range");
  return neighbors_[siteIndex];
}

Vector3i Cube::GetRelativePosition(size_t siteIndex) const
{
  if (siteIndex >= relativeIndexMatrix_.cols())
    throw std::out_of_range("Site index out of range");

  Vector3i pos = relativeIndexMatrix_.col(siteIndex);

  return pos;
}

Matrix3i Cube::GetBasis() const
{
  return basis_;
}

Matrix3Xi Cube::GetRelativePositionMatrix() const
{
  return relativeIndexMatrix_;
}

size_t Cube::GetSizeOfCube() const
{
  return sizeOfCube_;
}

size_t Cube::GetNumOfSites() const
{
  return sizeOfCube_ * sizeOfCube_ * sizeOfCube_;
}

size_t Cube::GetCentralSiteId() const
{
  int ci = sizeOfCube_ / 2;
  int cj = sizeOfCube_ / 2;
  int ck = sizeOfCube_ / 2;

  return ci * sizeOfCube_ * sizeOfCube_ + cj * sizeOfCube_ + ck;
}
