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

// From Atom basis class

RowVectorXd AtomBasis::GetAtomBasis(const Element &elementType)
{
  // Chebyshev basis
  if (basisType_ == "Chebyshev")
  {
    int index = -1;
    int i = 0;
    for (auto element : elementSet_)
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

    if (elementSet_.size() == 1)
      xVal = 0;
    else
      xVal = -1.0 + 2.0 * index / (elementSet_.size() - 1.0);

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
  else if (basisType_ == "Occupation")
  {
    RowVectorXd basis(elementSet_.size());
    basis.setZero();
    int i = 0;
    for (auto element : elementSet_)
    {
      if (element == elementType)
      {
        basis(i) = 1;
        return basis;
      }
      i++;
    }

    throw std::invalid_argument(
        "GetAtomBasis Error: Element '" + elementType.GetElementString() + "' not found in the element set.");
  }

  else
  {
    throw std::invalid_argument("GetAtomBasis Error: Invalid basis, use 'Chebyshev' or 'Occupation'.");
  }
}

RowVectorXd AtomBasis::GetCachedAtomBasis(const Element &elementType)
{
  RowVectorXd atomBasis;
  const string key = elementType.GetElementString();

  if (atomBasisHashMap_.Get(key, atomBasis))
  {

    // cout << "From cache" << endl;
    return atomBasis;
  }

  atomBasis = GetAtomBasis(elementType);

  atomBasisHashMap_.Add(key, atomBasis);

  // cout << "Added to the cache" << endl;

  return atomBasis;
}

RowVectorXd AtomBasis::GetCachedTensorProduct(
    const string &elementCluster,
    const vector<RowVectorXd> &basisVector,
    bool isSymmetric)
{
  RowVectorXd tensorProduct;
  if (tensorProductHashMap_.Get(elementCluster, tensorProduct))
  {
    // cout << "From cache" << endl;
    return tensorProduct;
  }

  tensorProduct = GetTensorProduct(basisVector, isSymmetric);

  tensorProductHashMap_.Add(elementCluster, tensorProduct);

  // cout << "Added to cache" << endl;

  return tensorProduct;
}

RowVectorXd AtomBasis::GetTensorProduct(const vector<RowVectorXd> &basisVector,
                                        bool isSymmetric)
{
  if (basisVector.empty())
    return RowVectorXd();

  // If only one basis vector, return it directly
  if (basisVector.size() == 1)
    return basisVector[0];

  // For symmetric product: average over all permutations
  // Cluster which are formed by sites which are equivalent due to symmetry
  // W-Ta and Ta-W will be equivalent if they are formed by sites which are
  // equivalent.
  //
  // W-Ta : [0 0 1 0]
  // Ta-W : [0 1 0 0]
  //
  // Due to symmetry basis vector for W-Ta will be
  // 1/2 * (basis(W)*basis(Ta) + basis(Ta)*basis(W))
  // which will result into [0 0.5 0.5 0]
  if (isSymmetric)
  {
    vector<int> indices(basisVector.size());
    for (size_t i = 0; i < indices.size(); ++i)
      indices[i] = i;

    std::sort(indices.begin(), indices.end());

    RowVectorXd sum;
    int count = 0;

    do
    {
      vector<RowVectorXd> permuted;
      for (int idx : indices)
        permuted.push_back(basisVector[idx]);

      RowVectorXd prod = GetTensorProduct(permuted, false);

      if (sum.size() == 0)
        sum = RowVectorXd::Zero(prod.size());

      sum += prod;
      ++count;

    } while (std::next_permutation(indices.begin(), indices.end()));

    return sum / static_cast<double>(count);
  }

  // Standard (non-symmetric) tensor product
  // Cluster which are formed by sites which are not equivalent.
  // W-Ta and Ta-W will not be equivalent
  // W-Ta : [0 0 1 0]
  // Ta-W : [0 1 0 1]
  RowVectorXd result = basisVector[0];

  for (size_t i = 1; i < basisVector.size(); ++i)
  {
    const RowVectorXd &vec = basisVector[i];
    RowVectorXd temp(result.size() * vec.size());

    for (int j = 0; j < result.size(); ++j)
    {
      temp.segment(j * vec.size(), vec.size()) = result(j) * vec;
    }

    result = temp;
  }

  return result;
}