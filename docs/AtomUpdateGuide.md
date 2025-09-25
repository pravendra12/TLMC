# UpdateAtomVector

## Description
Update the atom vector and atomic index vector of the `TiledSupercell`.

This function takes a vector of atomic indices (`atomicIndicesVector`) and updates
the internal atom representation (`atomVector_`) as well as the atomic index vector
(`atomIndexVector_`) of the supercell. Each element in `atomicIndicesVector`
corresponds to the atomic number of an element at a particular site in the
supercell, following the same ordering as the `TiledSupercell`'s internal atom vector.

### Important Notes
- The size of `atomicIndicesVector` must exactly match the total number of sites
  (`totalNumOfSites_`) in the supercell.
- The ordering of elements in `atomicIndicesVector` must correspond to the order
  used by this `TiledSupercell` object when it was dumped using
  `WriteAtomicIndicesToFile()` or a similar function.
  Using a vector from another supercell or in a different order may result in an
  inconsistent atom configuration.
- Each atomic index should represent a valid atomic number that can be used to
  construct an `Element` object. Passing invalid numbers may throw exceptions
  depending on the implementation of the `Element` class.

## Parameters
| Parameter | Type | Description |
|-----------|------|-------------|
| `atomicIndicesVector` | `std::vector<uint64_t>` | Vector of atomic indices representing the atomic numbers for each site in the supercell. |

## Exceptions
- Throws `std::runtime_error` if the size of the vector does not match the supercell.

## Example
```cpp
std::vector<uint64_t> atomicIndices = TiledSupercell::ReadAtomicIndicesFromFile("atom_indices.bin.gz");
supercell.UpdateAtomVector(atomicIndices);
