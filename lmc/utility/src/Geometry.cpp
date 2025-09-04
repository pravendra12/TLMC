
#include "Geometry.h" // Assumes Structure class is defined here

// Convert Cartesian to Fractional positions
Eigen::MatrixXd get_scaled_positions(
    const Eigen::MatrixXd &positions,
    const Eigen::Matrix3d &cell,
    bool wrap = true,
    const std::array<bool, 3> &pbc = {true, true, true})
{
  Eigen::MatrixXd fractional = cell.transpose().inverse() * positions.transpose();
  fractional.transposeInPlace();

  if (wrap)
  {
    for (int i = 0; i < 3; ++i)
    {
      if (pbc[i])
      {
        for (int j = 0; j < fractional.rows(); ++j)
        {
          // Apply modulo twice to handle negative values correctly
          fractional(j, i) = std::fmod(std::fmod(fractional(j, i), 1.0) + 1.0, 1.0);
        }
      }
    }
  }
  return fractional;
}

// Convert Fractional to Cartesian positions
Eigen::MatrixXd fractional_to_cartesian(
    const Structure &structure,
    const Eigen::MatrixXd &frac_positions)
{
  return frac_positions * structure.getCell();
}

// Convert Structure to Spglib cell
void structure_to_spglib_cell(
    const Structure &structure,
    double lattice[3][3],
    std::vector<std::array<double, 3>> &positions,
    std::vector<int> &numbers)
{
  // Copy cell
  Eigen::Matrix3d cell = structure.getCell();
  if (cell.determinant() < 1e-10)
  {
    throw std::runtime_error("Invalid lattice: determinant too small.");
  }
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      lattice[i][j] = cell(i, j);
    }
  }

  // Resize vectors
  const size_t natoms = structure.size();
  positions.resize(natoms);
  numbers.resize(natoms);

  // Get fractional positions
  Eigen::MatrixXd cart_positions(natoms, 3);
  for (size_t i = 0; i < natoms; ++i)
  {
    cart_positions.row(i) = structure.positionByIndex(i);
  }

  std::vector<bool> pbc_vec = structure.getPBC();
  std::array<bool, 3> pbc_array;

  for (size_t i = 0; i < 3; ++i)
    pbc_array[i] = pbc_vec[i];

  // Now you can call your function
  Eigen::MatrixXd frac_positions = get_scaled_positions(cart_positions, cell, true, pbc_array);

  // Fill positions and numbers
  auto atomic_numbers = structure.getAtomicNumbers();
  auto ptr = atomic_numbers.data();
  for (size_t i = 0; i < natoms; ++i)
  {
    positions[i] = {frac_positions(i, 0), frac_positions(i, 1), frac_positions(i, 2)};
    numbers[i] = ptr[i];
  }
}

// Get primitive structure using Spglib
Structure get_primitive_structure(
    const Structure &structure,
    bool no_idealize = true,
    bool to_primitive = true,
    double symprec = 1e-5)
{
  // Copy structure data
  double lattice[3][3];
  std::vector<std::array<double, 3>> positions;
  std::vector<int> numbers;
  structure_to_spglib_cell(structure, lattice, positions, numbers);
  int natoms = static_cast<int>(structure.size());

  // Standardize to primitive cell
  int n_atoms_primitive = spg_standardize_cell(
      lattice,
      reinterpret_cast<double (*)[3]>(positions.data()),
      numbers.data(),
      natoms,
      to_primitive ? 1 : 0,
      no_idealize ? 1 : 0,
      symprec);

  if (n_atoms_primitive == 0)
  {
    throw std::runtime_error("Spglib failed to find primitive cell, maybe reduce symprec.");
  }

  // Resize vectors
  positions.resize(n_atoms_primitive);
  numbers.resize(n_atoms_primitive);

  // Copy lattice
  Eigen::Matrix3d prim_cell;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      prim_cell(i, j) = lattice[i][j];
    }
  }

  // Copy positions (fractional)
  Eigen::MatrixXd frac_positionList(n_atoms_primitive, 3);
  for (int i = 0; i < n_atoms_primitive; ++i)
  {
    frac_positionList(i, 0) = std::round(positions[i][0] * 1e12) / 1e12;
    frac_positionList(i, 1) = std::round(positions[i][1] * 1e12) / 1e12;
    frac_positionList(i, 2) = std::round(positions[i][2] * 1e12) / 1e12;
  }

  // Convert to Cartesian
  Eigen::MatrixXd positionList = fractional_to_cartesian(structure, frac_positionList);
  pybind11::array_t<int> atomicNumbers = pybind11::array_t<int>(numbers.size());
  std::memcpy(atomicNumbers.mutable_data(), numbers.data(), n_atoms_primitive * sizeof(int));

  // Build primitive Structure
  Structure primitive(positionList, atomicNumbers, prim_cell, structure.getPBC());

  return primitive;
}

// and the corresponding chemical symbols for each site
std::pair<Structure, std::vector<std::vector<std::string>>> get_occupied_primitive_structure(
    const Structure &structure,
    const std::vector<std::vector<std::string>> &allowed_species,
    double symprec = 1e-5)
{
  // Determine if allowed_species is a flat list
  bool is_flat_list = !allowed_species.empty() &&
                      std::all_of(allowed_species.begin(), allowed_species.end(),
                                  [](const auto &inner)
                                  { return inner.size() == 1; }) &&
                      allowed_species.size() != structure.size();

  // Expand allowed_species to match structure size if flat list
  std::vector<std::vector<std::string>> expanded_species;
  if (is_flat_list)
  {
    std::vector<std::string> flat_species;
    for (const auto &inner : allowed_species)
    {
      flat_species.push_back(inner[0]);
    }
    expanded_species.resize(structure.size(), flat_species);
  }
  else
  {
    if (structure.size() != allowed_species.size())
    {
      throw std::runtime_error("Structure and allowed_species must have the same size unless allowed_species is a flat list.");
    }
    expanded_species = allowed_species;
  }

  // Get unique sorted symbols (sublattices)
  std::set<std::vector<std::string>> sorted_symbols_set;
  for (const auto &sym : expanded_species)
  {
    std::vector<std::string> sorted_sym = sym;
    std::sort(sorted_sym.begin(), sorted_sym.end());
    sorted_symbols_set.insert(sorted_sym);
  }
  std::vector<std::vector<std::string>> sorted_symbols(sorted_symbols_set.begin(), sorted_symbols_set.end());

  // Map sublattices to temporary atomic numbers using Element(chemicalSymbol).GetAtomicIndex()
  std::map<std::vector<std::string>, int> sublattice_to_atomic_number;
  for (const auto &sorted_sym : sorted_symbols)
  {
    if (sorted_sym.empty())
    {
      throw std::runtime_error("Empty symbol list in sublattice.");
    }
    // Use the first symbol in the sorted list to assign a unique atomic number
    Element element(sorted_sym[0]);
    int atomic_number = element.GetAtomicIndex();
    sublattice_to_atomic_number[sorted_sym] = atomic_number;
  }

  // Create a copy of the structure and assign temporary atomic numbers
  Structure decorated_primitive = structure;
  std::vector<int> new_numbers(structure.size());
  for (size_t i = 0; i < expanded_species.size(); ++i)
  {
    std::vector<std::string> sorted_sym = expanded_species[i];
    std::sort(sorted_sym.begin(), sorted_sym.end());
    auto it = sublattice_to_atomic_number.find(sorted_sym);
    if (it == sublattice_to_atomic_number.end())
    {
      throw std::runtime_error("Sublattice not found in mapping.");
    }
    new_numbers[i] = it->second; // Assign atomic number based on sublattice
  }
  decorated_primitive.setAtomicNumbers(py::cast(new_numbers));

  // Get primitive structure
  decorated_primitive = get_primitive_structure(decorated_primitive, true, true, symprec);

  // Map atomic numbers back to chemical symbols
  std::vector<std::vector<std::string>> primitive_chemical_symbols;
  auto prim_numbers = decorated_primitive.getAtomicNumbers();
  auto prim_numbers_ptr = prim_numbers.data();
  for (size_t i = 0; i < decorated_primitive.size(); ++i)
  {
    int atomic_number = prim_numbers_ptr[i];
    // Find the sublattice corresponding to this atomic number
    auto it = std::find_if(sublattice_to_atomic_number.begin(), sublattice_to_atomic_number.end(),
                           [atomic_number](const auto &pair)
                           { return pair.second == atomic_number; });
    if (it == sublattice_to_atomic_number.end())
    {
      throw std::runtime_error("Atomic number in primitive structure not found in sublattice mapping.");
    }
    primitive_chemical_symbols.push_back(it->first);
  }

  // Restore original expanded_species order where possible
  for (const auto &symbols : expanded_species)
  {
    std::vector<std::string> sorted_sym = symbols;
    std::sort(sorted_sym.begin(), sorted_sym.end());
    auto it = std::find_if(primitive_chemical_symbols.begin(), primitive_chemical_symbols.end(),
                           [&sorted_sym](const auto &s)
                           {
                             std::vector<std::string> sorted_s = s;
                             std::sort(sorted_s.begin(), sorted_s.end());
                             return sorted_s == sorted_sym;
                           });
    if (it != primitive_chemical_symbols.end())
    {
      *it = symbols; // Replace with original (unsorted) symbols
    }
  }

  /*
  // Wrap positions in the primitive structure
  Eigen::MatrixXd prim_positions(decorated_primitive.size(), 3);
  for (size_t i = 0; i < decorated_primitive.size(); ++i)
  {
    prim_positions.row(i) = decorated_primitive.positionByIndex(i);
  }
  Eigen::MatrixXd wrapped_positions = get_scaled_positions(
      prim_positions, decorated_primitive.getCell(), true, decorated_primitive.getPBC());
  decorated_primitive.setPositions(fractional_to_cartesian(decorated_primitive, wrapped_positions));
  */

  return {decorated_primitive, primitive_chemical_symbols};
}

/*
std::vector<std::string> get_wyckoff_sites(
    Structure structure,
    const std::vector<std::vector<std::string>> &map_occupations = {},
    double symprec = 1e-5,
    bool include_representative_atom_index = false)
{
  size_t natoms = structure.size();

  // Handle map_occupations
  std::vector<int> numbers(natoms);
  auto original_numbers = structure.getAtomicNumbers();
  auto ptr = original_numbers.data();

  for (size_t i = 0; i < natoms; ++i)
    numbers[i] = ptr[i];

  if (!map_occupations.empty())
  {
    if (map_occupations.size() == 1 && map_occupations[0].empty())
    {
      std::fill(numbers.begin(), numbers.end(), 1); // all H
    }
    else
    {
      for (size_t i = 0; i < natoms; ++i)
      {
        std::string symb = Element(numbers[i]).GetElementString();
        bool found = false;
        for (const auto &group : map_occupations)
        {
          if (std::find(group.begin(), group.end(), symb) != group.end())
          {
            numbers[i] = Element("H").GetAtomicIndex();
            found = true;
            break;
          }
        }
        // if not found, keep original
      }
    }
    structure.setAtomicNumbers(numbers);
  }

  // Convert to Spglib cell
  double lattice[3][3];
  std::vector<std::array<double, 3>> positions;
  std::vector<int> spg_numbers;
  structure_to_spglib_cell(structure, lattice, positions, spg_numbers);

  // Get symmetry dataset
  SpglibDataset dataset;
  int ret = spg_get_symmetry_dataset(lattice,
                                     reinterpret_cast<double (*)[3]>(positions.data()),
                                     spg_numbers.data(),
                                     natoms,
                                     symprec,
                                     &dataset);
  if (ret == 0)
    throw std::runtime_error("Spglib failed to get symmetry dataset.");

  // Number of unit cells
  Eigen::Matrix3d trans_mat;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      trans_mat(i, j) = dataset.transformation_matrix[i][j];
  double n_unitcells = std::abs(trans_mat.determinant());

  // Compute Wyckoff positions
  std::map<int, std::string> wyckoffs;
  std::vector<int> equivalent_atoms(dataset.n_atoms);
  for (int i = 0; i < dataset.n_atoms; i++)
    equivalent_atoms[i] = dataset.equivalent_atoms[i];

  for (int index : std::set<int>(equivalent_atoms.begin(), equivalent_atoms.end()))
  {
    int multiplicity = std::count(equivalent_atoms.begin(), equivalent_atoms.end(), index);
    multiplicity = static_cast<int>(std::round(multiplicity / n_unitcells));
    std::string wyckoff = std::to_string(multiplicity) + dataset.wyckoffs[index];
    if (include_representative_atom_index)
      wyckoff += "-" + std::to_string(index);
    wyckoffs[index] = wyckoff;
  }

  // Map Wyckoff positions to atoms
  std::vector<std::string> result(natoms);
  for (size_t i = 0; i < natoms; i++)
    result[i] = wyckoffs[equivalent_atoms[i]];

  spg_free_dataset(&dataset);
  return result;
}
*/