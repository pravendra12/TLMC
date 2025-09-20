#pragma once
#include <cstddef>

struct LatticeSiteEncodedMapping
{
  size_t latticeId;            // ID of the lattice in the supercell
  int encodedSmallConfigId; // encoded representation of the small config
};