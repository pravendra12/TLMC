#include <iostream>
#include <vector>
#include <stdexcept>
#include <tuple>
#include <Eigen/Dense>
#include <spglib.h>

using namespace std;
using namespace Eigen;

class MatrixOfEquivalentPositions {
private:
    MatrixXd translations;
    vector<Matrix3d> rotations;
    vector<vector<Vector3d>> positions;
    size_t n_symmetries;

public:
    MatrixOfEquivalentPositions(const MatrixXd& translations_input,
                                 const vector<Matrix3d>& rotations_input)
    {
        if (translations_input.rows() != (int)rotations_input.size()) {
            throw invalid_argument("The number of translations must equal the number of rotations.");
        }
        n_symmetries = rotations_input.size();
        translations = translations_input;
        rotations = rotations_input;
    }

    void build(const MatrixXd& fractional_positions) {
        size_t n_atoms = fractional_positions.rows();
        positions.resize(n_atoms);

        for (size_t i = 0; i < n_atoms; ++i) {
            positions[i].resize(n_symmetries);
            Vector3d atom_pos = fractional_positions.row(i);
            for (size_t j = 0; j < n_symmetries; ++j) {
                Vector3d new_pos = rotations[j] * atom_pos + translations.row(j).transpose();
                positions[i][j] = new_pos;
            }
        }
    }

    const vector<vector<Vector3d>>& get_equivalent_positions() const {
        return positions;
    }
};


#include <Eigen/Dense>
#include <spglib.h>
#include <tuple>
#include <vector>
#include <memory>
#include <stdexcept>

std::tuple<Eigen::MatrixXd, std::vector<Eigen::Matrix3d>> get_symmetry_from_spglib(
    const Eigen::Matrix3d &lattice,
    const Eigen::MatrixXd &positions,
    const std::vector<int> &atomic_numbers,
    double symprec)
{
    // Validate inputs
    if (positions.rows() != static_cast<Eigen::Index>(atomic_numbers.size()))
    {
        throw std::runtime_error("Number of positions must match number of atomic types.");
    }

    int num_atoms = positions.rows();

    // Convert lattice to C array
    double spg_lattice[3][3];
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            spg_lattice[i][j] = lattice(i, j);
        }
    }

    // Allocate positions in C format
    double (*spg_positions)[3] = new double[num_atoms][3];
    for (int i = 0; i < num_atoms; ++i)
    {
        spg_positions[i][0] = positions(i, 0);
        spg_positions[i][1] = positions(i, 1);
        spg_positions[i][2] = positions(i, 2);
    }

    // Copy atomic types
    std::vector<int> spg_types = atomic_numbers;

    // Get symmetry dataset
    SpglibDataset *dataset = spg_get_dataset(spg_lattice, spg_positions, spg_types.data(), num_atoms, symprec);
    delete[] spg_positions; // Free positions after use

    if (!dataset)
    {
        throw std::runtime_error("Failed to get symmetry dataset from spglib.");
    }

    // RAII cleanup for dataset
    struct DatasetDeleter
    {
        void operator()(SpglibDataset *ptr) const { spg_free_dataset(ptr); }
    };
    std::unique_ptr<SpglibDataset, DatasetDeleter> dataset_guard(dataset);

    // Extract rotations and translations
    std::vector<Eigen::Matrix3d> rotations(dataset->n_operations);
    Eigen::MatrixXd translations(dataset->n_operations, 3);

    for (int i = 0; i < dataset->n_operations; ++i)
    {
        Eigen::Matrix3d rot;
        for (int r = 0; r < 3; ++r)
        {
            for (int c = 0; c < 3; ++c)
            {
                rot(r, c) = static_cast<double>(dataset->rotations[i][r][c]);
            }
        }
        rotations[i] = rot;

        translations.row(i) << dataset->translations[i][0],
                               dataset->translations[i][1],
                               dataset->translations[i][2];
    }

    return std::make_tuple(translations, rotations);
}

// Main function equivalent to Python: matrix_of_equivalent_positions_from_structure
tuple<MatrixOfEquivalentPositions, MatrixXd>
matrix_of_equivalent_positions_from_structure(const Matrix3d& lattice,
                                              const MatrixXd& fractional_positions,
                                              const vector<int>& atomic_numbers,
                                              double symprec)
{
    auto [translations, rotations] = get_symmetry_from_spglib(lattice, fractional_positions, atomic_numbers, symprec);

    MatrixOfEquivalentPositions matrix_equiv(translations, rotations);
    matrix_equiv.build(fractional_positions);

    return make_tuple(matrix_equiv, fractional_positions);
}


