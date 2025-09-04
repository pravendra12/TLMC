#include "TestingFunctions.h"

void TestEnergyChangePredictor(Config &config, const string &predictorFilename)
{
  cout << "\n-------------------------------------------------------\n"
       << endl;
  cout << "Test Function for dE Computation for a cyclic path" << endl;
  EnergyPredictor energyPredictor(predictorFilename, config);

  vector<size_t> cyclicPath;
  bool pathFound = false;

  // Helper lambda: tries to find a valid cyclic path starting from startId
  auto findCyclicPathFromSite = [&config](size_t startId, vector<size_t> &outPath) -> bool
  {
    auto firstNN = config.GetNeighborLatticeIdVectorOfLattice(startId, 1);
    for (auto firstId : firstNN)
    {
      auto secondNN = config.GetNeighborLatticeIdVectorOfLattice(firstId, 1);
      for (auto secondId : secondNN)
      {
        if (secondId == startId)
          continue;

        auto thirdNN = config.GetNeighborLatticeIdVectorOfLattice(secondId, 1);
        for (auto thirdId : thirdNN)
        {
          auto fourthNN = config.GetNeighborLatticeIdVectorOfLattice(thirdId, 1);
          for (auto fourthId : fourthNN)
          {
            if (fourthId != startId)
              continue;

            vector<size_t> candidatePath = {startId, firstId, secondId, thirdId, fourthId};

            // Reject if intermediate elements not same or equal to start element
            auto intermediateRef = config.GetElementOfLattice(candidatePath[1]);
            bool intermediateSame = true;
            for (size_t i = 1; i < candidatePath.size() - 1; ++i)
            {
              if (config.GetElementOfLattice(candidatePath[i]) != intermediateRef)
              {
                intermediateSame = false;
                break;
              }
            }

            if (intermediateSame && intermediateRef != config.GetElementOfLattice(candidatePath[0]))
            {
              if (firstId != thirdId) // optional extra condition
              {
                outPath = candidatePath;
                return true;
              }
            }
          }
        }
      }
    }
    return false;
  };

  // Iterate over all sites until a valid cyclic path is found
  for (size_t siteId = 0; siteId < config.GetNumAtoms(); ++siteId)
  {
    if (findCyclicPathFromSite(siteId, cyclicPath))
    {
      pathFound = true;
      break;
    }
  }

  if (!pathFound)
  {
    cout << "No valid cyclic path found in the entire configuration." << endl;
    return;
  }

  // ------------------------------
  // Print the cyclic path nicely
  // ------------------------------
  cout << "Chosen cyclic path: ";
  print1DVector(cyclicPath);

  auto printElementOnPath = [](Config &config, const vector<size_t> &cyclicPath)
  {
    for (size_t i = 0; i < cyclicPath.size(); i++)
    {
      size_t latticeId = cyclicPath[i];
      auto elem = config.GetElementOfLattice(latticeId).GetElementString();
      cout << elem << "(" << latticeId << ")";
      if (i != cyclicPath.size() - 1)
        cout << " -> ";
    }
    cout << endl;
  };

  cout << "Elements along path: ";
  printElementOnPath(config, cyclicPath);

  // ------------------------------
  // Compute incremental energy changes along the cyclic path
  // ------------------------------
  cout << "\nIncremental energy changes along the cyclic path:" << endl;

  double totaldE = 0;
  for (size_t i = 0; i < cyclicPath.size() - 1; ++i)
  {
    size_t fromId = cyclicPath[i];
    size_t toId = cyclicPath[i + 1];

    cout << "Hop " << i << ": " << fromId << " -> " << toId << endl;

    cout << "Before Swap: ";
    printElementOnPath(config, cyclicPath);

    double energyBefore = energyPredictor.ComputeEnergyOfSite(config, fromId) +
                          energyPredictor.ComputeEnergyOfSite(config, toId);

    // Perform the hop (swap elements)
    config.LatticeJump(make_pair(fromId, toId));

    cout << "After  Swap: ";
    printElementOnPath(config, cyclicPath);

    double energyAfter = energyPredictor.ComputeEnergyOfSite(config, fromId) +
                         energyPredictor.ComputeEnergyOfSite(config, toId);

    double dE = energyAfter - energyBefore;
    cout << "ΔE = " << dE << endl;
    cout << "--------------------" << endl;

    totaldE += dE;
  }

  cout << "Total Energy Change: " << totaldE << endl;
  cout << "\n-------------------------------------------------------\n" << endl;
}
