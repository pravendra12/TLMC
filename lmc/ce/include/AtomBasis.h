/**
 * @file AtomBasis.h
 * @brief Provides functionality to compute basis vectors for atoms using different basis types.
 *
 * This header defines the `GetAtomBasis` function, which computes a basis vector
 * for a given element type based on the specified basis type ("Chebyshev" or "Occupation").
 * The function supports two types of basis:
 * - **Chebyshev Basis**: Computes Chebyshev polynomials of the first kind for the given element.
 * - **Occupation Basis**: Represents the presence of an element in a set as a one-hot encoded vector.
 */


#ifndef LMC_CE_INCLUDE_ATOMBASIS_H_
#define LMC_CE_INCLUDE_ATOMBASIS_H_

#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <unordered_map>
#include <unordered_set>
#include "Element.hpp"

using namespace std;
using namespace Eigen;

/**
 * @brief Generates a basis vector for a specified element and basis type.
 *
 * @param elementType The element for which the basis vector is generated.
 * @param elementSet A set of elements providing the context for basis computation.
 * @param basisType The type of basis to generate. Supported options:
 *                  - "Chebyshev": Generates Chebyshev polynomials of the first kind.
 *                  - "Occupation": Generates a one-hot encoded vector.
 * @return A RowVectorXd representing the generated basis vector.
 *
 * @throws std::invalid_argument If the element is not part of the element set.
 * @throws std::invalid_argument If the basis type is invalid.
 *
 * Chebyshev Basis:
 * - Computes Chebyshev polynomials of the first kind, T_n(x), where `x` is derived
 *   from the element's position in the element set.
 * - The value of `x` is scaled to the range [-1, 1] based on the element's index
 *   in the set and the total number of elements.
 * - The first polynomial T_0(x) is always 1, and T_1(x) is equal to `x`.
 * - Higher-order polynomials are computed using the recurrence relation:
 *   T_n(x) = 2 * x * T_(n-1)(x) - T_(n-2)(x).
 * - Currently, the maximum polynomial degree is set to 1.
 *
 * Occupation Basis:
 * - Produces a one-hot encoded vector with a 1 at the position corresponding to
 *   the specified element, and 0s elsewhere.
 */
RowVectorXd GetAtomBasis(const Element &elementType,
                         const set<Element> &elementSet,
                         const string &basisType);


#endif // LMC_CE_INCLUDE_ATOMBASIS_H_