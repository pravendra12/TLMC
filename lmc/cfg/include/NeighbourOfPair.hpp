#pragma once
#include <cstddef>
#include "LatticeSiteEncodedMapping.hpp"

// Represents a neighbor lattice site relative to a pair of lattice sites.
// This struct is used when working with lattice pairs, where each neighbor
// must be associated with one of the two lattice sites in the pair.
// It also stores the encoded configuration index for decoding purposes.

struct NeighbourOfPair
{
    LatticeSiteEncodedMapping latticeSiteInfo; // Mapping information for the neighboring lattice site,
                                               // including its lattice ID and encoded configuration index.

    // Enum indicating which lattice in the pair this neighbor is associated with.
    // Using a strongly typed enum ensures type safety and prevents invalid assignments.
    // Possible values:
    //   First  = neighbor belongs to the first lattice in the pair
    //   Second = neighbor belongs to the second lattice in the pair
    enum class SourceSite : uint8_t
    {
        First = 0,
        Second = 1
    };

    // Stores the source lattice in the pair that this neighbor belongs to.
    // Using the strongly typed enum prevents accidental invalid assignments
    // and makes the code more readable.
    SourceSite origin;
};
