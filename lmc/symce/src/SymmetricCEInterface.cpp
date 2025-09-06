#include "SymmetricCEInterface.h"

SymmetricCEInterface::SymmetricCEInterface(
    const Config &config,
    const vector<string> &allowedElements,
    const vector<double> &clusterCutoffs,
    const double symprec) : pyGuard_{}
{
  // Start Python Interpreter
  // py::scoped_interpreter guard{};

  try
  {
    // Python Utilities
    fs::path pyPath = fs::path(__FILE__).parent_path().parent_path() / "py";

    py::module_ sys = py::module_::import("sys");
    sys.attr("path").attr("append")(pyPath.string());

    // Import Python modules
    pyStandalone_ = py::module_::import("standalone");
    py::module_ pyASEAtoms = py::module_::import("ase.atoms");

    // ASE Atoms
    pyAtomsObject_ = ConvertConfigToASEAtoms(
        pyASEAtoms,
        config);

    py::print(pyAtomsObject_);

    // Cutoffs
    for (const auto &cutoff : clusterCutoffs)
    {
      pyCutoffs_.append(cutoff);
    }

    vector<int> atomicNumberVector;

    for (const auto &element : allowedElements)
    {
      pyAllowedSpecies_.append(py::str(element));

      atomicNumberVector.emplace_back(
        Element(element).GetAtomicIndex()
      );
    }

    allowedAtomicNumbersVector_ = {atomicNumberVector};

    symprec_ = symprec;
  }
  catch (const py::error_already_set &e)
  {
    cerr << "Python error: " << e.what() << endl;
    throw runtime_error("Failed to initialize Python module");
  }
  catch (const exception &e)
  {
    cerr << "C++ error: " << e.what() << endl;
    throw;
  }
}

py::object SymmetricCEInterface::ConvertConfigToASEAtoms(
    py::module_ &aseAtoms,
    const Config &config)
{
  // Convert Config to Structure
  auto structure = ConvertConfigToStructure(config);

  // Positions : Convert Eigen matrices to NumPy arrays
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> pos = structure.getPositions();
  py::array_t<double> positions({pos.rows(), pos.cols()}, pos.data());

  // Basis
  Eigen::Matrix3d basis = structure.getCell();
  py::array_t<double> cell({3, 3}, basis.data());

  // Atomic Numbers
  py::array_t<int> atomicNumbers = structure.getAtomicNumbers();

  // PBC
  py::list pbc;
  for (bool b : structure.getPBC())
    pbc.append(b);

  // Create ASE Atoms object
  py::object Atoms = aseAtoms.attr("Atoms");
  py::object atomObject = Atoms(
      py::arg("positions") = positions,
      py::arg("numbers") = atomicNumbers,
      py::arg("cell") = cell,
      py::arg("pbc") = pbc);

  return atomObject;
}

Structure SymmetricCEInterface::GetPrimitiveStructure()
{
  auto pyPrimitiveStructure = pyStandalone_.attr("GetInputsForOrbitListPrimStructre")(
      pyAtomsObject_,
      pyAllowedSpecies_,
      pyCutoffs_,
      symprec_);

  cout << "Primitive Structure" << endl;
  py::print(pyPrimitiveStructure);

  // Basis
  py::array_t<double> pyCell = pyPrimitiveStructure.attr("get_cell")();
  auto cellBuffer = pyCell.unchecked<2>(); // 2D array

  Eigen::Matrix3d cell;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      cell(i, j) = cellBuffer(i, j);

  // Positions
  py::array_t<double> pyPositions = pyPrimitiveStructure.attr("get_positions")();
  auto positionBuffer = pyPositions.unchecked<2>();

  size_t numAtoms = pyPositions.shape(0);
  Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor> positions(numAtoms, 3);

  for (size_t i = 0; i < numAtoms; i++)
  {
    for (size_t j = 0; j < 3; j++)
    {
      positions(i, j) = positionBuffer(i, j);
    }
  }

  // PBC
  py::list pyPbc = pyPrimitiveStructure.attr("get_pbc")();
  vector<bool> pbc(3);
  for (int i = 0; i < 3; i++)
    pbc[i] = pyPbc[i].cast<bool>();

  // Atomic Numbers List
  py::array_t<int> pyAtomicNumbers = pyPrimitiveStructure.attr("get_atomic_numbers")();

  Structure primitiveStructure(
      positions,
      pyAtomicNumbers,
      cell,
      pbc);

  primitiveStructure.setAllowedAtomicNumbers(allowedAtomicNumbersVector_);

  return primitiveStructure;
}

vector<vector<LatticeSite>> SymmetricCEInterface::GetMatrixOfEquivalentSites()
{
  auto pyMatrixOfEquivalentSites = pyStandalone_.attr("GetInputsForOrbitListPMLatticeSites")(
      pyAtomsObject_,
      pyAllowedSpecies_,
      pyCutoffs_,
      symprec_);

  vector<vector<LatticeSite>> matrixOfEquivalentSites;

  for (auto rowObject : pyMatrixOfEquivalentSites)
  {
    vector<LatticeSite> rowVector;

    for (auto item : rowObject)
    {
      auto tupleItem = item.cast<py::tuple>(); // (int, tuple(int,int,int))
      int id = tupleItem[0].cast<int>();

      auto vectorTuple = tupleItem[1].cast<py::tuple>();
      Eigen::Vector3i vec3(vectorTuple[0].cast<int>(),
                           vectorTuple[1].cast<int>(),
                           vectorTuple[2].cast<int>());

      rowVector.emplace_back(id, vec3);
    }
    matrixOfEquivalentSites.push_back(move(rowVector));
  }

  return matrixOfEquivalentSites;
}

vector<vector<vector<LatticeSite>>> SymmetricCEInterface::GetNeighbourLists()
{
  auto pyNeighborLists = pyStandalone_.attr("GetInputsForOrbitListNeighbourList")(
      pyAtomsObject_,
      pyAllowedSpecies_,
      pyCutoffs_,
      symprec_);

  vector<vector<vector<LatticeSite>>> neighborLists;

  for (auto rowObject : pyNeighborLists)
  {
    py::list row = rowObject.cast<py::list>();
    vector<vector<LatticeSite>> rowVector;
    for (auto colObject : row)
    {
      py::list col = colObject.cast<py::list>();
      vector<LatticeSite> colVector;
      for (auto item : col)
      {
        auto tupleItem = item.cast<py::tuple>(); // (int, tuple(int,int,int))
        int id = tupleItem[0].cast<int>();

        auto vectorTuple = tupleItem[1].cast<py::tuple>();
        Eigen::Vector3i vec3(vectorTuple[0].cast<int>(),
                             vectorTuple[1].cast<int>(),
                             vectorTuple[2].cast<int>());

        colVector.emplace_back(id, vec3);
      }
      rowVector.push_back(colVector);
    }
    neighborLists.push_back(rowVector);
  }

  return neighborLists;
}

double SymmetricCEInterface::GetPositionTolerance()
{
  auto pyPositionTolerance = pyStandalone_.attr("GetInputsForOrbitListTolerance")(
      pyAtomsObject_,
      pyAllowedSpecies_,
      pyCutoffs_,
      symprec_);

  return pyPositionTolerance.cast<double>();
}

double SymmetricCEInterface::GetFractionalPositionTolerance()
{
  auto pyFractionalPositionTolerance = pyStandalone_.attr("GetInputsForOrbitListFractionalTolerance")(
      pyAtomsObject_,
      pyAllowedSpecies_,
      pyCutoffs_,
      symprec_);

  return pyFractionalPositionTolerance.cast<double>();
}
