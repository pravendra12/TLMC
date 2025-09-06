#include "InterfaceUtility.h"

Config ConvertStructureToConfig(const Structure &structure)
{
  // Basis
  Matrix3d basis = structure.getCell();

  MatrixXd positions = structure.getPositions();         // Nx3
  Matrix3Xd relativePositionMatrix(3, positions.rows()); // 3 x N

  // Relative Positions
  for (Index i = 0; i < positions.rows(); ++i)
  {
    Vector3d pos = positions.row(i);
    relativePositionMatrix.col(i) = pos;
  }

  // Atomic Numbers
  auto atomicNumbers = structure.getAtomicNumbers(); // py::array_t<int>
  auto buf = atomicNumbers.request();                // buffer info
  int *ptr = static_cast<int *>(buf.ptr);
  size_t numLattices = buf.shape[0];

  std::vector<Element> atomVector(numLattices);
  for (size_t i = 0; i < numLattices; ++i)
  {
    // Declaring element vector using atomic number
    atomVector[i] = Element(ptr[i]); 
  }

  // Build Config
  return Config(basis, relativePositionMatrix, atomVector);
}

Structure ConvertConfigToStructure(const Config &config)
{
  size_t numLattices = config.GetNumLattices();

  // Positions
  Matrix<double, Dynamic, 3, RowMajor> positions(numLattices, 3);
  vector<int> atomicNumberVector(numLattices);

  for (size_t i = 0; i < numLattices; ++i)
  {
    auto pos = config.GetCartesianPositionOfLattice(i);
    positions.row(i) = pos;
    atomicNumberVector[i] = config.GetElementOfLattice(i).GetAtomicIndex();
  }

  py::array_t<int> atomicNumbersPy(static_cast<py::ssize_t>(numLattices));
  auto buf = atomicNumbersPy.mutable_data();
  copy(atomicNumberVector.begin(), atomicNumberVector.end(), buf);

  // Basis or cell
  Matrix3d cell = config.GetBasis();

  // PBC : All simulations are mostly under pbc
  vector<bool> pbc(3, true); 

  // Construct Structure
  Structure structure(positions, atomicNumbersPy, cell, pbc);
  return structure;
}
