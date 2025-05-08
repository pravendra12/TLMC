#include "AtomBasis.h"

RowVectorXd GetAtomBasis(const Element &elementType,
                         const set<Element> &elementSet,
                         const string &basisType)
{
  // Chebyshev basis
  if (basisType == "Chebyshev")
  {
    int index = -1;
    int i = 0;
    for (auto element : elementSet)
    {
      if (element == elementType)
      {
        index = i;
        break;
      }
      i++;
    }

    if (index == -1)
    {
      throw std::invalid_argument("Element not in element set");
    }

    double xVal;

    if (elementSet.size() == 1)
      xVal = 0;
    else
      xVal = -1.0 + 2.0 * index / (elementSet.size() - 1.0);

    int maxDegree = 1;

    int size = maxDegree + 1;

    RowVectorXd basis(size);
    basis(0) = 1.0; // T_0(x)

    if (maxDegree >= 1)
    {
      basis(1) = xVal; // T_1(x)
      for (int n = 2; n <= maxDegree; ++n)
      {
        basis(n) = 2 * xVal * basis(n - 1) - basis(n - 2); // T_n(x)
      }
    }

    return basis;
  }

  // Occupancy basis
  else if (basisType == "Occupation")
  {
    RowVectorXd basis(elementSet.size());
    basis.setZero();
    int i = 0;
    for (auto element : elementSet)
    {
      if (element == elementType)
      {
        basis(i) = 1;
        return basis;
      }
      i++;
    }

    throw std::invalid_argument("Element not in element set");
  }

  else
  {
    throw std::invalid_argument("Error: Invalid basis, use 'Chebyshev' or 'Occupation'.");
  }
}