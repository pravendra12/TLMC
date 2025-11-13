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

#include "GenerateClosedLoops.h"

int main()
{
  auto cfg = Config::GenerateSupercell(
      10, 3.2, "Mo", "BCC");

  cfg.UpdateNeighborList({3});

  // first neighbour
  auto neighbourList = cfg.GetNeighborLists()[0];

  auto closedLoops = GenerateClosedLoops(
      neighbourList,
      50,
      50);

  std::string baseOutputDir = "/home/pravendra3/Documents/TLMC/bin/testingClosedLoops";

  for (size_t i = 0; i < closedLoops.size(); ++i)
  {
    const auto &loop = closedLoops[i];

    // Print the loop
    std::cout << "Loop " << i << " (length " << loop.size() << "): ";
    for (size_t j = 0; j < loop.size(); ++j)
    {
      std::cout << loop[j];
      if (j + 1 < loop.size())
        std::cout << " -> ";
    }
    std::cout << std::endl;

    // Write step-wise configs
    WriteLoopStepConfigs(cfg, loop, baseOutputDir, i);

    std::cout << "Wrote loop " << i << " with " << loop.size() << " steps in folder: "
              << baseOutputDir << "/loop_" << i << std::endl
              << std::endl;
  }
}
