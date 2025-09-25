#ifndef LMC_CFG_INCLUDE_TILEDSUPERCELL_H_
#define LMC_CFG_INCLUDE_TILEDSUPERCELL_H_

#include "Cube.h"
#include "Config.h"
#include "Element.hpp"
#include "Constants.hpp"
#include "Eigen/Dense"
#include <iostream>
#include <iomanip>
#include "LatticeSiteMapping.hpp"
#include "LatticeSiteEncodedMapping.hpp"
#include "NeighbourOfPair.hpp"
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/filesystem.hpp>
#include <vector>
#include <fstream>
#include <cstdint>

using namespace std;
using namespace Eigen;
namespace io = boost::iostreams;

// This class represents the superconfig formed by placing the
// config inside cube madeup of those config

class TiledSupercell
{
public:
  // Constructor
  TiledSupercell(
      Config &smallConfig,
      const Cube &cubeObj);

  // Get the small config or prim which was used to build the tiledsupercell
  const Config &GetSmallConfig() const;

  const Matrix3d &GetSuperBasis() const;

  // Get the cube object which was used to build teh tiledSupercell
  const Cube &GetCube() const;

  size_t GetNumOfSitesPerSmallConfig() const;
  size_t GetNumOfSmallConfig() const;
  size_t GetTotalNumOfSites() const;

  // Returns the neighbour list
  // <LATTICE_ID, ENCODED_SMALL_CFG_INDEX>
  // returns teh neighbuor <nnLatticeId, encodedSmallConfigId>
  // encodedSmallConfigId can be used to get the correct index using the sorted
  // neighobur of the current smallConfig in which latticeId lies
  vector<LatticeSiteEncodedMapping> GetNeighborLatticeIdVectorOfLattice(
      const size_t &latticeId,
      const size_t &distanceOrder) const;

  // return neighbours upto max bond order
  vector<LatticeSiteEncodedMapping> GetNeighborLatticeIdsUpToOrder(
      const size_t &latticeId,
      const size_t &maxBondOrder) const;

  // It returns the neighbouring LatticeSiteEncodedMapping for a latticeIdPair
  // Now it is independent of the configIdx for the two sites because the
  // neighbouring sites contain the information in encoded form regarding which
  // supercell they belong to
  // For more details check the comments in the function
  vector<NeighbourOfPair> GetNeighboringLatticeIdSetOfPair(
      const pair<size_t, size_t> &latticeIdPair,
      const size_t &maxBondOrder,
      const bool removeLatticeIdDuplicates) const;

  // these latticePair are bascially teh
  // <LatticeId, smallConfigId> which will be used to
  // get the atomId then the swap will happen
  void LatticeJump(
      const pair<LatticeSiteMapping, LatticeSiteMapping> &latticeSiteMappingPair);

  // Sets element of a site give the latticeAndConfigIdPair
  // latticeAndConfigIdPair : <LATTICE_ID, SMALL_CONFIG_IDX>
  // SMALL_CONFIG_IDX is the index of config inside the cube
  void SetElementAtSite(
      const LatticeSiteMapping &latticeSiteMapping,
      const Element &element);

  // Returns the element at site give the latticeAndConfigIdPair
  // latticeAndConfigIdPair : <LATTICE_ID, SMALL_CONFIG_IDX>
  // SMALL_CONFIG_IDX is the index of config inside the cube [NOT THE ENCODED_INDEX]
  Element GetElementAtSite(
      const LatticeSiteMapping &latticeSiteMapping) const;

  LatticeSiteMapping GetVacancySiteId() const;

  /**
   * @brief Decodes an encoded lattice site mapping into an actual lattice site mapping.
   *
   * Each encoded neighbor stores:
   *   - `latticeId` : ID of the lattice site within a small configuration.
   *   - `encodedSmallConfigId` : Index of the neighboring small configuration
   *                              relative to the reference small configuration.
   *                              A value of -1 indicates the neighbor is in the same
   *                              small configuration as the reference.
   *
   * This function converts the encoded representation to a full `LatticeSiteMapping`
   * containing the actual lattice ID and small configuration ID.
   *
   * @param latticeSiteEncoding Encoded lattice site information.
   * @param refLatticeSiteMapping The reference lattice site mapping used to resolve
   *                              the encoded small configuration index.
   * @return LatticeSiteMapping The corresponding lattice site mapping with actual IDs.
   */
  LatticeSiteMapping GetLatticeSiteMappingFromEncoding(
      const LatticeSiteEncodedMapping &latticeSiteEncoding,
      const LatticeSiteMapping &refLatticeSiteMapping) const;

  Vector3d GetRelativePositionOfLatticeSiteMapping(
      const LatticeSiteMapping &latticeSiteId) const;

  Vector3d GetRelativeDistanceVectorLattice(
      const LatticeSiteMapping &site1,
      const LatticeSiteMapping &site2) const;

  /**
   * @brief Convert a global atom ID to its corresponding smallConfig index and
   * lattice ID.
   *
   * This function maps a single atom's global index in the large tiled
   * supercell to the specific smallConfig it belongs to and the lattice site
   * within that config.
   *
   * @param atomId  The global atom index in the supercell.
   * @return        A pair {latticeId, smallConfigId} where:
   * LATTICE_ID : Lattice index of the atom inside the smallConfig.
   * SMALL_CONFIG_IDX : Index of the smallConfig in the tiled supercell.
   * @throws std::out_of_range if atomId exceeds the total number of atoms in
   * the supercell.
   */
  LatticeSiteMapping GetLatticeSiteMappingFromAtomId(
      const size_t &atomId) const;

  /**
   * @brief Convert a smallConfig index and lattice ID to the global atom ID.
   *
   * This function computes the global atom index in the large supercell given
   * the index of the smallConfig and the lattice site within that config.
   *
   * @param latticeAndConfigIdPair  <LATTICE_ID, SMALL_CONFIG_IDX>
   * LATTICE_ID : Lattice index of the atom inside the smallConfig.
   * SMALL_CONFIG_IDX : Index of the smallConfig in the tiled supercell.
   * @return               Global atom index in the supercell.
   * @throws std::out_of_range if the computed atom index exceeds the total number of atoms.
   */
  size_t GetAtomIdFromLatticeAndConfigId(
      const LatticeSiteMapping &latticeSiteMapping) const;

  // Builds a map from latticeId -> neighbors for all sites in the supercell
  // LATTICE_ID_IN_SMALL_CFG : {<NN_LATTICE_ID_IN_SMALL_CFG , SMALL_CFG_ID_IN_CUBE>}
  // User will provide a maxBondOrder such that the smallCfg must have its neighbours
  // being updated // need to add some check for the same
  void UpdateNeighbourLists(const size_t maxBondOrder);

  // Maps the atom information to a 1-D vector
  // this has been made public so that one can load a config or atom vector
  // then assign it to the tiledSuprecell
  void UpdateAtomVector(const Config &config);



  /**
   * @brief Update the TiledSupercell's atom vector using a dumped atomic indices vector.
   * @see docs/AtomUpdateGuide.md
   */
  void UpdateAtomVector(const vector<uint64_t> &atomicIndicesVector);

  // I/O

  Config MakeSupercell() const;

  static Config MakeSupercellFromAtomInfo(const Config &smallConfig,
                                          const Cube &cubeObj,
                                          const vector<Element> &atomVector);

  // Writes the atom info to a text file
  void WriteAtomVectorInfoToFile(const string &filename) const;

  // Writes the atom info to binary file
  void WriteAtomVectorInfoToBinary(const string &filename) const;

  // Write the Atom Index Vector to a binary compressed file
  void WriteAtomicIndicesToFile(const string &filename) const;

  // to read the atom vector from file
  static vector<Element> ReadAtomVectorInfoFromFile(const string &filename);

  static vector<Element> ReadAtomVectorInfoFromBinary(const string &filename);

  static vector<uint64_t> ReadAtomicIndicesFromFile(const string &filename);

private:
  const Config &smallConfig_; // Dont want smallConfig to be const so as to update the neighbours inside the TiledSupercell
  const Cube cubeObj_;

  const Matrix3d superBasis_{};

  vector<Element> atomVector_{};
  // This vector will be used for dumping the atomic information for efficieny purpose
  vector<uint64_t> atomIndexVector_{};

  // neighbourList_ is a 3-level nested container that stores neighbor information
  // for each lattice site in the tiled supercell.
  //
  // Structure:
  //   neighbourList_[bondOrder][siteId][k] = { neighborSiteId, neighborConfigIdx }
  //
  // - Level 1 (outer vector): groups neighbors by bond order
  //     * index 0 → first nearest neighbors (1nn)
  //     * index 1 → second nearest neighbors (2nn)
  //     * index 2 → third nearest neighbors (3nn)
  //     * ...
  //
  // - Level 2 (middle vector): for each site in the supercell, store its neighbors
  //     * neighbourList_[bondOrder][siteId] is a vector of all neighbors of
  //       the given site at that bond order.
  //
  // - Level 3 (inner vector of pairs): list of neighbors for the given site
  //     * Each entry is a pair { neighborSiteId, configIdx }
  //         - neighborSiteId : lattice index of the neighbor
  //         - configIdx      : which smallConfig (inside the cube) contains it
  //
  // Example:
  //   neighbourList_[0][5] = vector of all first nearest neighbors of site 5.
  //   neighbourList_[1][10][0] = { neighborSiteId, configIdx } for the first
  //                              second-nearest neighbor of site 10.
  vector<vector<vector<LatticeSiteEncodedMapping>>> neighbourList_{};

  const size_t numSitesPerSmallConfig_{};
  // Total smallConfig in the cube
  const size_t numSmallConfig_{};

  // Total sites : numSmallConfig_ * numSitesPerSmallConfig_
  const size_t totalNumOfSites_{};

  // Helper methods

  // This function checks whether relative to a reference lattice Id the
  // its neibhouringLatticeId whether it lies in same small config or it was
  // neighbour due to PBC
  bool LatticeInSameConfig(
      const size_t &referenceLatticeId,
      const size_t &neighbouringLatticeId) const;

  /**
   * @brief Map a neighbor lattice site to the correct smallConfig index inside the tiled cube.
   *
   * This function takes a lattice site (`referenceSiteId`) in a given smallConfig
   * and one of its neighbor lattice sites (`neighborSiteId`). If the neighbor
   * also belongs to the same smallConfig (without crossing periodic boundaries),
   * then the owning config index is trivial (`configIdxInCube` itself).
   *
   * However, in a tiled cube of smallConfigs, neighbors may fall into adjacent
   * smallConfigs due to periodic boundary conditions (PBC). In that case, we need
   * to determine which neighboring smallConfig inside the cube actually contains
   * the neighbor.
   *
   * The algorithm works as follows:
   *   1. Compute the scaled (fractional) position of the reference site
   *      by adding its cube offset and dividing by cube size.
   *   2. For each candidate neighboring config index returned by
   *      cubeObj.GetNeighbors(configIdxInCube):
   *        - Compute the scaled position of the neighbor site as if it belonged
   *          to that candidate config.
   *        - Apply the minimum image convention (wrap into [-0.5, 0.5]) to
   *          handle PBC in fractional coordinates.
   *        - Compute the Euclidean distance between the reference and neighbor.
   *   3. Select the config index that yields the minimum wrapped distance.
   *
   * @param smallConfig       The unit cell / smallConfig definition.
   * @param referenceSiteId   Lattice site index of the reference atom in smallConfig.
   * @param neighborSiteId    Lattice site index of the neighbor in smallConfig.
   * @param cubeObj           The Cube object describing how smallConfigs tile the supercell.
   * @param configIdxInCube   Index of the current smallConfig instance inside the cube.
   *
   * @return size_t           Index of the smallConfig inside the cube
   *                          that actually contains the neighbor site.
   */
  // Finds which neighboring smallConfig (inside the cube) contains the neighbor site
  // that is closest to the given reference site, considering periodic boundary conditions
  // for the cube.
  size_t FindNearestNeighborConfigIndex(
      const size_t referenceLatticeId,     // site in the current small config
      const size_t neighborLatticeId,      // site whose owning config index must be determined
      const size_t configIdxInCube) const; // index of the current small config in which reference Site lies

  // Returns a list of neighbors for a lattice site in a given small config within the tiled supercell.
  // Each pair: <neighborLatticeId, smallConfigIndexInCube>
  // <LATTICE_ID_IN_SMALL_CFG , SMALL_CFG_ID_IN_CUBE>
  vector<pair<size_t, size_t>> GetLatticeNeighbors(
      const size_t latticeId,       // Lattice site index in the small config
      const size_t configIdxInCube, // Index of this small config in the cube
      const size_t maxBondOrder     // Number of neighbor shells
  ) const;

  void InitializeAtomVector();

  void PrintTiledSupercell() const;
};

#endif // LMC_CFG_INCLUDE_TILEDSUPERCELL_H_