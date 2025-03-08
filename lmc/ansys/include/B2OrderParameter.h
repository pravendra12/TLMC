#ifndef LMC_ANSYS_INCLUDE_B2ORDERPARAMETER_H_
#define LMC_ANSYS_INCLUDE_B2ORDERPARAMETER_H_

#include "Config.h"
#include <omp.h>
#include <utility>
#include <unordered_set>
using namespace std;

/**
 * @class B2OrderParameter
 * @brief Computes the B2 order parameter for a given element in a configuration.
 * 
 * The B2 order parameter is a measure of atomic ordering in binary alloys, particularly in BCC structures.
 * It quantifies the difference in site occupancy of a given element between two distinct sublattices (alpha and beta).
 */
class B2OrderParameter
{
  public:
    /**
     * @brief Constructs a B2OrderParameter object and initializes lattice site classifications.
     * 
     * This constructor initializes the alpha and beta lattice sites by identifying their positions
     * in the given configuration. It calls `InitializeAlphaLatticeSites()` and `InitializeBetaLatticeSites()`
     * to classify the sites into two groups based on second nearest neighbor (NN) relationships.
     * 
     * @param config The configuration object containing lattice information, including element occupancy
     *               and neighbor relationships.
     */
    B2OrderParameter(const Config &config);


    /** 
     * @brief Computes the B2 order parameter for a given element.
     * 
     * The B2 order parameter measures the degree of ordering in a B2-type structure.
     * It is defined as:
     * 
     * \f[
     * B2 = \frac{f_{\alpha} - f_{\beta}}{f_{\alpha} + f_{\beta}}
     * \f]
     * 
     * where:
     *
     * - \f$f_{\alpha}\f$  is the fractional occupancy of the element at alpha sites.
     * - \f$f_{\alpha}\f$  is the fractional occupancy of the element at beta sites.
     * 
     * A perfectly ordered B2 structure has \( B2 = 1 \), while a completely disordered structure has \( B2 = 0 \).
     * 
     * @param element The element for which the B2 order parameter is computed.
     * @return The B2 order parameter value for the specified element.
     */  
    double GetB2OrderParameter(const Element &element);

    /**
     * @brief Computes the fractional occupancy of an element at alpha lattice sites.
     * 
     * The fractional occupancy is calculated as:
     * \f[
     * f_{\alpha} = \frac{\text{Number of alpha sites occupied by the element}}{\text{Total number of alpha sites}}
     * \f]
     * This function iterates through the alpha sites and counts how many are occupied by the given element.
     * 
     * @param element The element whose occupancy is being calculated.
     * @return Fractional occupancy of the element at alpha sites (a value between 0 and 1).
     */
    double GetAlphaSiteOccupancy(const Element &element);
    
    /**
     * @brief Computes the fractional occupancy of an element at beta lattice sites.
     * 
     * Similar to `GetAlphaSiteOccupancy()`, this function calculates the fraction of beta sites occupied
     * by the given element using the formula:
     * \f[
     * f_{\beta} = \frac{\text{Number of beta sites occupied by the element}}{\text{Total number of beta sites}}
     * \f]
     * It iterates through the beta sites and counts how many are occupied by the given element.
     * 
     * @param element The element whose occupancy is being calculated.
     * @return Fractional occupancy of the element at beta sites (a value between 0 and 1).
     */
    double GetBetaSiteOccupancy(const Element &element);

    /**
     * @brief Retrieves the set of lattice sites classified as alpha sites.
     * 
     * Alpha sites are a subset of the lattice that are identified based on second nearest
     * neighbor relationships. These sites are used to compute the B2 order parameter.
     * 
     * @return An unordered set containing the indices of alpha lattice sites.
     */
    unordered_set<size_t> GetAlphaLatticeSites();
    
    /**
     * @brief Retrieves the set of lattice sites classified as beta sites.
     * 
     * Beta sites are defined as the complement of alpha sites in the lattice.
     * These sites are used in the calculation of the B2 order parameter.
     * 
     * @return An unordered set containing the indices of beta lattice sites.
     */
    unordered_set<size_t> GetBetaLatticeSites();
        
  private:
    
    /**
     * @brief Identifies and initializes lattice sites as alpha sites based on second nearest neighbor (NN) relationships.
     * 
     * The function classifies alpha sites by analyzing the second nearest neighbor list.
     * It follows these steps:
     * 1. Retrieves the second nearest neighbor list from the configuration.
     * 2. Iterates over all sites, checking their second nearest neighbors.
     * 3. Ensures that each new alpha site does not have a first nearest neighbor (NN1) 
     *    bond with an already classified alpha site, maintaining the required sublattice separation.
     * 4. Stores the identified alpha sites in an unordered set.
     * 
     * This method effectively partitions the lattice into two interpenetrating sublattices,
     * which are essential for measuring B2 ordering.
     */
    void InitializeAlphaLatticeSites();

    /**
     * @brief Identifies and initializes lattice sites as beta sites, complementing the alpha sites.
     * 
     * This function assigns all remaining lattice sites (not classified as alpha sites)
     * to the beta sublattice. Since beta sites are simply the complement of alpha sites
     * in the lattice, this function:
     * 1. Iterates through all lattice site indices.
     * 2. Checks if a site is already classified as an alpha site.
     * 3. If not, assigns it as a beta site.
     * 
     * This ensures that every lattice site is categorized as either an alpha or a beta site,
     * enabling the computation of the B2 order parameter.
     */
    void InitializeBetaLatticeSites();
    
    /**
     * @brief Reference to the configuration object containing lattice information.
     * 
     * This object provides details about lattice sites, neighbor relationships, and element occupancy.
     * It is used to classify sites into alpha and beta sublattices and to compute the B2 order parameter.
     */
    const Config &config;

    /**
     * @brief Total number of lattice sites in the configuration.
     * 
     * Retrieved from the configuration using `config.GetNumLattices()`, this variable stores
     * the total count of lattice sites, which is useful for iterating over all sites
     * during initialization and analysis.
     */
    const size_t numLattice = config.GetNumLattices();
    
    /**
     * @brief Set of lattice sites classified as alpha sites.
     * 
     * Alpha sites are one of the two sublattices in a B2-ordered structure.
     * These sites are determined based on second nearest neighbor (NN2) relationships.
     * They are used to compute the fractional occupancy for calculating the B2 order parameter.
     */
    unordered_set<size_t> alphaLatticeSites;

    /**
     * @brief Set of lattice sites classified as beta sites.
     * 
     * Beta sites form the complementary sublattice to alpha sites in a B2-ordered structure.
     * Once alpha sites are identified, the remaining lattice sites are assigned to the beta sublattice.
     * The fraction of an element occupying beta sites contributes to the B2 order parameter calculation.
     */
    unordered_set<size_t> betaLatticeSites;

};






#endif //LMC_ANSYS_INCLUDE_B2ORDERING_H_