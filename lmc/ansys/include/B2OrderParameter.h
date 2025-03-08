#ifndef LMC_ANSYS_INCLUDE_B2ORDERPARAMETER_H_
#define LMC_ANSYS_INCLUDE_B2ORDERPARAMETER_H_

#include "Config.h"
#include <omp.h>
#include <utility>
#include <unordered_set>
using namespace std;

class B2OrderParameter
{
  public:
    // Constructor
    B2OrderParameter(const Config &config);

    double GetB2OrderParameter(const Element &element);

    double GetAlphaSiteOccupancy(const Element &element);

    double GetBetaSiteOccupancy(const Element &element);

    // Get Alpha Sites
    unordered_set<size_t> GetAlphaLatticeSites();
    
    // Get Beta Sites
    unordered_set<size_t> GetBetaLatticeSites();
    
    
  private:
    
    // Initializes alpha lattice sites
    void InitializeAlphaLatticeSites();

    // Initializes beta lattice sites
    void InitializeBetaLatticeSites();
        
    const Config &config;
    const size_t numLattice = config.GetNumLattices();
    
    // Alpha Sites
    unordered_set<size_t> alphaLatticeSites;

    // Beta Sites
    unordered_set<size_t> betaLatticeSites;

};






#endif //LMC_ANSYS_INCLUDE_B2ORDERING_H_