#ifndef LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_
#define LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_

#include <omp.h>
#include <iostream>
#include <stdexcept>
#include <filesystem>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "SubLatticeOccupancy.h"
#include "TiledSupercell.h"
#include "ShortRangeOrderTLMC.h"
#include "B2ClusterTLMC.h"

using namespace std;
namespace fs = filesystem;

namespace ansys
{

  class Traverse
  {
  public:
    /**
     * @brief Construct a Traverse object
     *
     * @param initialSteps
     * @param incrementSteps
     * @param cutoffs
     * @param logType
     */
    Traverse(unsigned long long int initialSteps,
             unsigned long long int incrementSteps,
             const string logType);

    virtual ~Traverse();

    void RunAnsys(
        const TiledSupercell &tiledSupercell,
        const SubLatticeOccupancy &subLatticeOccupancy,
        const set<Element> &elementSet,
        const unordered_set<size_t> &convertToConfigSet) const;

    void RunAnsysOnConfig(
        const TiledSupercell &tiledSupercell,
        const vector<size_t> &atomicIndicesVector,
        const SubLatticeOccupancy &subLatticeOccupancy,
        const set<Element> &elementSet,
        ostringstream &oss,
        const bool &saveConfig = false,
        const string &outfilename = "") const;

  private:
    string GetHeaderFrameString(const set<Element> &elementSet) const;

    const unsigned long long initialSteps_;
    const unsigned long long incrementSteps_;
    unsigned long long finalSteps_;
    const string logType_;

    using MapVariant =
        variant<unordered_map<unsigned long long, double>, unordered_map<unsigned long long, string>>;

    unordered_map<string, MapVariant> logMap_;

    mutable ofstream frameOfs_;
  };

} // namespace ansys

#endif // LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_