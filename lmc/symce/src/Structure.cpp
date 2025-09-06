#include "Structure.hpp"

using namespace Eigen;

/**
  @details Initializes an icet Structure instance.
  @param positions list of positions in Cartesian coordinates
  @param atomicNumbers list of atomic numbers
  @param cell cell metric
  @param pbc periodic boundary conditions
**/
Structure::Structure(const Matrix<double, Dynamic, 3, RowMajor> &positions,
                     const py::array_t<int> &atomicNumbers,
                     const Matrix3d &cell,
                     const std::vector<bool> &pbc)
    : _atomicNumbers(atomicNumbers), _cell(cell), _pbc(pbc)
{
    _scaledPositions = positions * _cell.inverse();
}

/**
  @details Returns the position of a lattice site.
  @param latticeSite site for which to obtain the position
  @returns a 3-dimensional position vector
*/
Vector3d Structure::getPosition(const LatticeSite &latticeSite) const
{
    if (latticeSite.index() >= (size_t)_scaledPositions.rows())
    {
        std::ostringstream msg;
        msg << "Site index out of bounds";
        msg << " index: " << latticeSite.index();
        msg << " number of sites: " << _scaledPositions.rows();
        msg << " (Structure::getPosition)";
        throw std::out_of_range(msg.str());
    }
    Vector3d position = (_scaledPositions.row(latticeSite.index()) + latticeSite.unitcellOffset().transpose().cast<double>()) * _cell;
    return position;
}

/**
  @details Returns the position of a specific site in Cartesian coordinates.
  @param index index of the site
*/
Vector3d Structure::positionByIndex(const size_t &index) const
{
    if (index >= (size_t)_scaledPositions.rows())
    {
        std::ostringstream msg;
        msg << "Index out of bounds";
        msg << " index: " << index;
        msg << " number of sites: " << _scaledPositions.rows();
        msg << " (Structure::positionByIndex)";
        throw std::out_of_range(msg.str());
    }
    return _scaledPositions.row(index) * _cell;
}

/**
  @brief Returns the positions (in Cartesian coordinates) of all atoms in this structure.
*/
Matrix<double, Dynamic, 3, RowMajor> Structure::getPositions() const
{
    return _scaledPositions * _cell;
}

/**
  @details This function returns the LatticeSite object the position of
  which matches the input position to the tolerance specified for this
  structure.

  The algorithm commences by extracting the fractional position.
  From the fractional position the unitcell offset is taken by rounding the
  fractional coordinates to the nearest integer.
  By subtracting the unitcell offset from the fractional position and taking
  the dot product with the cell the position relative to the primitive cell is
  found.
  The index is found by searching for the remainder position in structure.
  If no index is found a runtime_error is thrown.

  @param position position to match in Cartesian coordinates
  @param fractionalPositionTolerance tolerance applied when comparing positions in fractional coordinates

  @returns LatticeSite object
*/
LatticeSite Structure::findLatticeSiteByPosition(const Vector3d &position, const double fractionalPositionTolerance) const
{
    Vector3d scaledPosition = _cell.transpose().partialPivLu().solve(position);
    /// Loop over all positions
    for (size_t i = 0; i < (size_t)_scaledPositions.rows(); i++)
    {
        Vector3d fractionalDistanceVector = scaledPosition - _scaledPositions.row(i).transpose();

        // Check whether whether the fractionalDistanceVector is all integers,
        // if it is we have found the corresponding lattice site
        Vector3d latticeVector = {round(fractionalDistanceVector[0]),
                                  round(fractionalDistanceVector[1]),
                                  round(fractionalDistanceVector[2])};
        if ((fractionalDistanceVector - latticeVector).norm() < fractionalPositionTolerance)
        {
            return LatticeSite(i, latticeVector.cast<int>());
        }
    }

    std::ostringstream msg;
    msg << "Failed to find site by position (findLatticeSiteByPosition)." << std::endl;
    msg << "Try increasing symprec or position_tolerance." << std::endl;
    msg << "position: " << position[0] << " " << position[1] << " " << position[2] << std::endl;
    msg << "scaled position: " << scaledPosition[0] << " " << scaledPosition[1] << " " << scaledPosition[2] << std::endl;
    msg << "fractional position tolerance: " << fractionalPositionTolerance;
    throw std::runtime_error(msg.str());
}

/**
  @brief Prescribes the atomic numbers that are allowed on each site in the structure.
  @param atomicNumbers
    Nested list with the atomic numbers allowed on each site in the structure.
**/
void Structure::setAllowedAtomicNumbers(const std::vector<std::vector<int>> &atomicNumbers)
{
    if (atomicNumbers.size() != size())
    {
        std::ostringstream msg;
        msg << "Size of input list incompatible with structure";
        msg << " length: " << atomicNumbers.size();
        msg << " nsites: " << size();
        msg << " (Structure::setAllowedAtomicNumbers)";
        throw std::out_of_range(msg.str());
    }
    _allowedAtomicNumbers = atomicNumbers;
}

/**
  @brief Returns the allowed atomic numbers per site in this structure.
*/
const std::vector<std::vector<int>> &Structure::allowedAtomicNumbers() const
{
    if (!hasAllowedAtomicNumbers())
    {
        std::ostringstream msg;
        msg << "Allowed atomic numbers per site not set in this structure";
        msg << " (Structure::allowedAtomicNumbers)";
        throw std::out_of_range(msg.str());
    }
    return _allowedAtomicNumbers;
}
