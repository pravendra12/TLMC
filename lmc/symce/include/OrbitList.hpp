#pragma once

#include <vector>

#include "Cluster.hpp"
#include "LatticeSite.hpp"
#include "ManyBodyNeighborList.hpp"
#include "Orbit.hpp"
#include "Structure.hpp"
#include "Symmetry.hpp"
#include "VectorOperations.hpp"
#include "VectorOperations.hpp"

/**
This class serves as a container for a sorted list of orbits and provides associated functionality.
*/

class OrbitList
{
public:
    /// Empty constructor.
    OrbitList(){};

    /// Constructor with structure only.
    OrbitList(const Structure &structure) { _structure = structure; };

    /// Constructs orbit list from a set of neighbor lists, a matrix of equivalent sites, and a structure.
    OrbitList(const Structure &,
              const std::vector<std::vector<LatticeSite>> &,
              const std::vector<std::vector<std::vector<LatticeSite>>> &,
              const double);

    /// Sort orbit list.
    void sort(const double);

    /// Adds an orbit.
    void addOrbit(const Orbit &orbit);

    /// Adding-to-existing (+=) operator.
    OrbitList &operator+=(const OrbitList &);

    /// Returns the orbit of the given index.
    const Orbit &getOrbit(unsigned int) const;

    // @todo Add description.
    void addColumnsFromMatrixOfEquivalentSites(std::vector<std::vector<std::vector<LatticeSite>>> &,
                                               std::unordered_set<std::vector<int>, VectorHash> &,
                                               const std::vector<int> &) const;

    /// Returns the number of orbits.
    size_t size() const { return _orbits.size(); }

    // Returns the first column of the matrix of equivalent sites.
    std::vector<LatticeSite> getReferenceLatticeSites() const;

    // Returns rows of the matrix of equivalent sites that match the lattice sites.
    std::vector<int> getReferenceLatticeSiteIndices(const std::vector<LatticeSite> &,
                                                    bool = true) const;

    // Returns true if the cluster includes at least on site from the unit cell at the origin.
    bool validCluster(const std::vector<LatticeSite> &) const;

    // @todo Add description.
    std::vector<LatticeSite> getTranslatedSites(const std::vector<LatticeSite> &,
                                                const unsigned int) const;

    /// @todo Add description.
    std::vector<std::vector<LatticeSite>> getSitesTranslatedToUnitcell(const std::vector<LatticeSite> &,
                                                                       bool sort = true) const;

    /// @todo Add description.
    std::vector<std::pair<std::vector<LatticeSite>, std::vector<int>>> getMatchesInMatrixOfEquivalentSites(const std::vector<std::vector<LatticeSite>> &) const;

    /// @todo Add description.
    void transformSiteToSupercell(LatticeSite &,
                                  const Structure &,
                                  std::unordered_map<LatticeSite, LatticeSite> &,
                                  const double) const;

    /// Returns symmetry related site groups via _matrixOfEquivalentSites and translations
    std::vector<std::vector<LatticeSite>> getSymmetryRelatedSiteGroups(const std::vector<LatticeSite> &) const;

    /// Returns the matrix of equivalent sites used to construct the orbit list.
    std::vector<std::vector<LatticeSite>> getMatrixOfEquivalentSites() const { return _matrixOfEquivalentSites; }

    /// Removes an orbit identified by its index.
    void removeOrbit(const size_t);

    /// Removes all orbits that have inactive sites.
    void removeOrbitsWithInactiveSites();

    /// Returns the orbits in this orbit list.
    const std::vector<Orbit> &orbits() const { return _orbits; }

    /// Returns the structure.
    const Structure &structure() const;

    /// Merge two orbits.
    void mergeOrbits(int index1, int index2);

private:
    /// Contains all the orbits in the orbit list.
    std::vector<Orbit> _orbits;

    /// Lattice sites used as references, equivalent to the first column of the _matrixOfEquivalentSites
    std::vector<LatticeSite> _referenceLatticeSites;

    /// Matrix of equivalent sites.
    std::vector<std::vector<LatticeSite>> _matrixOfEquivalentSites;

    /// Structure for which orbit list was constructed.
    Structure _structure;

    /// Create an orbit by deducing the proper permutations
    Orbit createOrbit(const std::vector<std::vector<LatticeSite>> &);
};
