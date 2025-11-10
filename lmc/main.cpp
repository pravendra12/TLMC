/*! \file  main.cpp
 *  \brief File for the main function.
 */

#include "Home.h"
/*
int main(int argc, char *argv[])
{
  if (argc == 1)
  {
    std::cout << "No input parameter filename." << std::endl;
    return 1;
  }
  api::Parameter parameter(argc, argv);
  api::Print(parameter);
  api::Run(parameter);
}
  */

#include "RateCorrector.hpp"

int main()
{
  // pred::RateCorrector rateCorrector(0.001, "/home/pravendra3/Documents/TLMC/bin/vfe_output.txt");

  size_t smallConfigSize = 10;
  double latticeParam = 3.2;
  auto smallCfg = Config::GenerateSupercell(
      smallConfigSize,
      latticeParam,
      "Mo",
      "BCC");
  size_t cubeSize = 2;
  Cube cubeObj(cubeSize);

  TiledSupercell tiledSupercell(
      smallCfg,
      cubeObj);

  auto cfg = Config::GenerateAlloySupercell(
      20,
      3.2,
      "BCC",
      {"Mo", "Ta"},
      {50, 50},
      1);

  cfg.SetElementOfLattice(0, Element("X"));

  auto concMap = tiledSupercell.GetConcentrationMap();

  for (const auto entry : concMap)
  {
    cout << entry.first.GetElementString() << " : " << entry.second << endl;
  }

  tiledSupercell.UpdateAtomVector(cfg);
  concMap = tiledSupercell.GetConcentrationMap();

  auto atomIndicesVector = tiledSupercell.GetAtomIndexVector();

  for (const auto entry : concMap)
  {
    cout << entry.first.GetElementString() << " : " << entry.second << endl;
  }

  tiledSupercell.UpdateAtomVector(atomIndicesVector);
  concMap = tiledSupercell.GetConcentrationMap();

  for (const auto entry : concMap)
  {
    cout << entry.first.GetElementString() << " : " << entry.second << endl;
  }

  RateCorrector rateCorrector(concMap, "/home/pravendra3/Documents/TLMC/bin/vfe_output.txt");

  for (int i = 0; i < 4; ++i)
  {
    auto start = std::chrono::high_resolution_clock::now();

    double corr = rateCorrector.GetTimeCorrectionFactor(1000.0);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Call " << i + 1
              << ": Time correction factor = " << corr
              << ", elapsed time = " << elapsed.count() << " s"
              << std::endl;
  }

  return 0;
}
