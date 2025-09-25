#include "SaveCubeAsConfig.h"

void SaveCubeAsConfig(
    const string &filename,
    const Cube &cubeObj,
    const TiledSupercell &tiledSupercell)
{
  vector<Element> atomVector;

  for (int i = 0; i < cubeObj.GetNumOfSites(); i++)
  {
    atomVector.emplace_back(Element("X"));
  }

  // Get integer positions from Cube
  Matrix3Xi intPositions = cubeObj.GetRelativePositionMatrix();

  // Convert to double
  Matrix3Xd positions = intPositions.cast<double>();

  // Convert to relative (fractional) positions in [0,1]
  Matrix3Xd relativePositions = positions / double(cubeObj.GetSizeOfCube());

  // Now use relativePositions to construct Config
  Config cubeConfig(
      tiledSupercell.GetSuperBasis(),
      relativePositions, // 3D fractional coordinates [0,1]
      atomVector);

  // Write config
  Config::WriteConfig(filename, cubeConfig);
}
