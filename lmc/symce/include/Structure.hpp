#pragma once

#include <Eigen/Dense>
#include <pybind11/numpy.h>

#include "LatticeSite.hpp"

using namespace Eigen;
namespace py = pybind11;

/**
  @brief Class for storing a structure.
  @details This class stores the cell metric, positions, atomic numbers, and
  periodic boundary conditions that describe a structure. It also holds
  information pertaining to the components that are allowed on each site and
  provides functionality for computing distances between sites.
*/
class Structure
{
public:
    /// Default constructor.
    Structure(){};

    /// Overloaded constructor.
    Structure(const Matrix<double, Dynamic, 3, RowMajor> &,
              const py::array_t<int> &,
              const Matrix3d &,
              const std::vector<bool> &);

    /// Return the position of a site in Cartesian coordinates.
    Vector3d getPosition(const LatticeSite &) const;

    /// Return the position of specific site in the structure in Cartesian coordinates.
    Vector3d positionByIndex(const size_t &) const;

    /// Returns positions of all sites in the structure in Cartesian coordinates.
    Matrix<double, Dynamic, 3, RowMajor> getPositions() const;

    /// Return LatticeSite object that matches the given position.
    LatticeSite findLatticeSiteByPosition(const Vector3d &, const double) const;

    /// Returns the size of the structure, i.e., the number of sites.
    size_t size() const { return _atomicNumbers.size(); }

    /// Set atomic numbers.
    void setAtomicNumbers(const py::array_t<int> &atomicNumbers) { _atomicNumbers = atomicNumbers; }

    /// Returns atomic numbers.
    const py::array_t<int> &getAtomicNumbers() const { return _atomicNumbers; }

    /// Returns periodic boundary conditions.
    std::vector<bool> getPBC() const { return _pbc; }

    /// Returns the cell metric.
    Matrix<double, 3, 3> getCell() const { return _cell; }

    /// Set allowed components for each site by vector.
    void setAllowedAtomicNumbers(const std::vector<std::vector<int>> &);

    /// Checks whether allowed atomic numbers are set
    bool hasAllowedAtomicNumbers() const { return _allowedAtomicNumbers.size() == size(); }

    /// Get allowed components for each site.
    const std::vector<std::vector<int>> &allowedAtomicNumbers() const;

private:
    /// List of atomic numbers.
    py::array_t<int> _atomicNumbers;

    /// Positions of sites in Cartesian coordinates.
    Matrix<double, Dynamic, 3, RowMajor> _scaledPositions;

    /// Cell metric.
    Matrix3d _cell;

    /// Periodic boundary conditions.
    std::vector<bool> _pbc;

    /// Specifies the atomic numbers that are allowed on each site.
    std::vector<std::vector<int>> _allowedAtomicNumbers;
};
