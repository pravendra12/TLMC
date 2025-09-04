#ifndef LMC_CE_INCLUDE_BASISSET_H_
#define LMC_CE_INCLUDE_BASISSET_H_

#include <set>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <unordered_map>
#include <unordered_set>
#include "Element.hpp"
#include "LruCache.h"
#include "AtomClusterType.hpp"
#include <boost/functional/hash.hpp>

using namespace std;
using namespace Eigen;

class BasisSet
{
public:
  BasisSet(
      const set<Element> &elementSet,
      const string &basisType,
      const bool &debug = false, 
      const size_t cacheCapacity = 50);

  VectorXd GetCachedTensorProduct(
      const AtomClusterType &atomClusterType);

  VectorXd GetCachedAtomBasis(
      const Element &elementType);

  VectorXd GetAtomBasis(
      const Element &elementType);

  VectorXd GetTensorProduct(
      const AtomClusterType &atomClusterType);

  double GetBasisProduct(const AtomClusterType &atomClusterType);

  double GetCachedBasisProduct(const AtomClusterType &atomClusterType);

private:
  void TestBasisSet();

  VectorXd ComputeChebyshevBasis(
      size_t index,
      size_t setSize,
      size_t degree = 1);

  VectorXd ComputeOccupationBasis(
      size_t index,
      size_t setSize);

  const set<Element> elementSet_;
  const string basisType_;

  // Hash key will be used as key for LRU
  LruCache<size_t, VectorXd> atomBasisHashMap_;
  LruCache<size_t, VectorXd> tensorProductHashMap_;
  LruCache<size_t, double> basisProductHashMap_;

};

#endif // LMC_CE_INCLUDE_BASISSET_H_