#include "BasisSet.h"

BasisSet::BasisSet(
    const set<Element> &elementSet,
    const string &basisType,
    const bool &debug,
    const size_t cacheCapacity) : elementSet_(elementSet),
                                  basisType_(basisType),
                                  atomBasisHashMap_(cacheCapacity),
                                  tensorProductHashMap_(cacheCapacity), 
                                  basisProductHashMap_(cacheCapacity)
{
  if (!(basisType_ == "Occupation" || basisType_ == "Chebyshev"))
  {
    throw runtime_error("Only Occupation or Chebyshev basis is supported, got " + basisType_);
  }

  if (debug)
  {
    TestBasisSet();
  }
}

VectorXd BasisSet::GetCachedTensorProduct(const AtomClusterType &atomClusterType)
{
  VectorXd tensorProduct;

  auto atomClusterHashKey = boost::hash<AtomClusterType>()(atomClusterType);

  if (tensorProductHashMap_.Get(atomClusterHashKey, tensorProduct))
  {
    return tensorProduct;
  }

  tensorProduct = GetTensorProduct(atomClusterType);

  tensorProductHashMap_.Add(atomClusterHashKey, tensorProduct);
  return tensorProduct;
}

VectorXd BasisSet::ComputeChebyshevBasis(
    size_t index,
    size_t setSize,
    size_t degree)
{
  VectorXd basis(degree + 1);
  // Scale index to [-1, 1]
  double x = setSize > 1 ? (2.0 * index / (setSize - 1) - 1.0) : 0.0;

  // Chebyshev polynomials of the first kind
  basis(0) = 1.0; // T_0(x) = 1
  if (degree >= 1)
  {
    basis(1) = x; // T_1(x) = x
    for (size_t n = 2; n <= degree; ++n)
    {
      // T_n(x) = 2x * T_(n-1)(x) - T_(n-2)(x)
      basis(n) = 2.0 * x * basis(n - 1) - basis(n - 2);
    }
  }
  return basis;
}

VectorXd BasisSet::ComputeOccupationBasis(
    size_t index,
    size_t setSize)
{
  /*
  Occupation Basis
  [1, 1, 1] // A
  [0, 1, 0] // B
  [0, 0, 1] // C

  Ref: Constructing multicomponent cluster expansions with machine-learning and chemical embedding

  */
  VectorXd basis = VectorXd::Zero(setSize);

  if (index == 0)
  {
    // Reference element (A) → all ones
    basis.setOnes();
  }
  else
  {
    // Other elements → one-hot at their position
    basis(index) = 1.0;
  }

  return basis;
}

VectorXd BasisSet::GetAtomBasis(
    const Element &elementType)
{
  auto it = elementSet_.find(elementType);
  if (it == elementSet_.end())
  {
    throw invalid_argument("Element" + elementType.GetElementString() + "not found in element set.");
  }

  size_t index = distance(elementSet_.begin(), it);
  size_t setSize = elementSet_.size();

  if (basisType_ == "Chebyshev")
  {
    return ComputeChebyshevBasis(index, setSize);
  }
  else if (basisType_ == "Occupation")
  {
    return ComputeOccupationBasis(index, setSize);
  }
}

VectorXd BasisSet::GetCachedAtomBasis(const Element &elementType)
{
  VectorXd atomBasis;

  auto elementHashKey = boost::hash<Element>()(elementType);

  if (atomBasisHashMap_.Get(elementHashKey, atomBasis))
  {
    return atomBasis;
  }

  // Add it to cache
  atomBasis = GetAtomBasis(elementType);
  atomBasisHashMap_.Add(elementHashKey, atomBasis);

  return atomBasis;
}

VectorXd BasisSet::GetTensorProduct(const AtomClusterType &atomClusterType)
{
  auto elementVector = atomClusterType.GetElementVector();

  if (elementVector.empty())
    return VectorXd(); // handle empty cluster

  // Start with the first element’s basis
  VectorXd result = GetCachedAtomBasis(elementVector[0]);
  if (elementVector.size() == 1)
    return result;

  // Multiply iteratively
  for (size_t i = 1; i < elementVector.size(); ++i)
  {
    const VectorXd &vec = GetCachedAtomBasis(elementVector[i]);
    VectorXd temp(result.size() * vec.size());

    for (int j = 0; j < result.size(); ++j)
      temp.segment(j * vec.size(), vec.size()) = result(j) * vec;

    result.swap(temp);
  }

  return result;
}

double BasisSet::GetBasisProduct(const AtomClusterType &atomClusterType)
{
  auto elementVector = atomClusterType.GetElementVector();

  if (elementVector.empty())
    return 1.0;

  double basisProduct = 1.0;

  for (size_t i = 0; i < elementVector.size(); ++i)
  {
    const VectorXd &vec = GetCachedAtomBasis(elementVector[i]);

    // Use last element of the vector (non-constant Chebyshev term)
    basisProduct *= vec(vec.size() - 1);
  }

  return basisProduct;
}

double BasisSet::GetCachedBasisProduct(const AtomClusterType &atomClusterType)
{
  double basisProduct;

  auto atomClusterHashKey = boost::hash<AtomClusterType>()(atomClusterType);

  if (basisProductHashMap_.Get(atomClusterHashKey, basisProduct))
  {
    return basisProduct;
  }

  basisProduct = GetBasisProduct(atomClusterType);

  basisProductHashMap_.Add(atomClusterHashKey, basisProduct);
  return basisProduct;
}
void BasisSet::TestBasisSet()
{
  // Set fixed-point notation with 3 decimal places for vectors
  cout << fixed << setprecision(3);

  // Header
  cout << "====================================\n";
  cout << "         Testing BasisSet\n";
  cout << "====================================\n\n";

  // Print element set
  cout << "Elements: ";
  bool first = true;
  for (const auto &element : elementSet_)
  {
    if (!first)
      cout << ", ";
    cout << element.GetElementString();
    first = false;
  }
  cout << "\n\n";

  // Print atom basis vectors
  cout << "--- Atom Basis ---\n";
  for (const auto &element : elementSet_)
  {
    const Eigen::VectorXd basis = GetAtomBasis(element);
    const Eigen::VectorXd cachedBasis = GetCachedAtomBasis(element);
    cout << left << setw(4) << element.GetElementString()
         << ": " << basis.transpose()
         << " | Cached: " << cachedBasis.transpose() << "\n";
  }
  cout << "\n";

  // Print tensor products and basis products for atom pairs
  cout << "--- Tensor Product and Basis Product (Pairs) ---\n";
  for (auto it1 = elementSet_.begin(); it1 != elementSet_.end(); ++it1)
  {
    for (auto it2 = it1; it2 != elementSet_.end(); ++it2)
    {
      const AtomClusterType pair(*it1, *it2);
      const Eigen::VectorXd tensor = GetTensorProduct(pair);
      const Eigen::VectorXd cached = GetCachedTensorProduct(pair);
      double basisProduct = GetCachedBasisProduct(pair);

      cout << left << setw(8)
           << it1->GetElementString() + "-" + it2->GetElementString()
           << ": Tensor: " << tensor.transpose()
           << " | Cached Tensor: " << cached.transpose()
           << " | Cached Basis Product: " << basisProduct
           << "\n";
    }
  }
  cout << "\n";

  // Print tensor products and basis products for atom triplets
  cout << "--- Tensor Product and Basis Product (Triplets) ---\n";
  for (auto it1 = elementSet_.begin(); it1 != elementSet_.end(); ++it1)
  {
    for (auto it2 = it1; it2 != elementSet_.end(); ++it2)
    {
      for (auto it3 = it2; it3 != elementSet_.end(); ++it3)
      {
        const AtomClusterType triplet(*it1, *it2, *it3);
        const Eigen::VectorXd tensor = GetTensorProduct(triplet);
        const Eigen::VectorXd cached = GetCachedTensorProduct(triplet);
        double basisProduct = GetCachedBasisProduct(triplet);

        cout << left << setw(10)
             << it1->GetElementString() + "-" + it2->GetElementString() + "-" + it3->GetElementString()
             << ": Tensor: " << tensor.transpose()
             << " | Cached Tensor: " << cached.transpose()
             << " | Cached Basis Product: " << basisProduct
             << "\n";
      }
    }
  }
  cout << "\n";

  // Print cache sizes
  cout << "Cache Sizes:\n";
  cout << "- Atom Basis Cache: " << atomBasisHashMap_.GetSizeOfCacheList() << "\n";
  cout << "- Tensor Product Cache: " << tensorProductHashMap_.GetSizeOfCacheList() << "\n";
  cout << "- Basis Product Cache: " << basisProductHashMap_.GetSizeOfCacheList() << "\n";

  // Footer
  cout << "\n====================================\n";
  cout << "      End of BasisSet Test\n";
  cout << "====================================\n";
}
