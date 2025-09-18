// main for Tiled KMC

#include "Config.h"
#include "Cube.h"
#include "TiledSupercell.h"

int main()
{
  size_t smallConfigSize = 5;
  double latticeParam = 3.2;
  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");

  vector<double> cutoffs = {3};
  smallCfg.UpdateNeighborList(cutoffs);

  // Now to make a 50x50x50 supercell with 10x10x10 small suprecells
  // One would need to stack 5 of these small cfg in a cube
  // total 5x5x5 cube = 125 cubes
  // Total Atoms = 125 * 2000
  Cube cubeObj(5);

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

  cout << tiledSupercell.GetElementAtSite({0, 0}).GetElementString() << endl;

  tiledSupercell.SetElementAtSite({0, 0}, Element("X"));

  cout << tiledSupercell.GetElementAtSite({0, 0}).GetElementString() << endl;

  auto centralIdInConfig = smallCfg.GetCentralAtomLatticeId();
  auto centralIdInCube = cubeObj.GetCentralSiteId();

  tiledSupercell.SetElementAtSite({centralIdInConfig, centralIdInCube}, Element("X"));

  cout << tiledSupercell.GetElementAtSite({centralIdInConfig, centralIdInCube}).GetElementString() << endl;

  return 0;
}