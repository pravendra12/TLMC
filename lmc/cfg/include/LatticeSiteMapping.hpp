#pragma once
#include <cstddef>

struct LatticeSiteMapping
{
  size_t latticeId;     // lattice ID in the supercell
  size_t smallConfigId; // ID of the small config it belongs to

  // Equality operator
  bool operator==(const LatticeSiteMapping &other) const
  {
    return latticeId == other.latticeId && smallConfigId == other.smallConfigId;
  }
};
