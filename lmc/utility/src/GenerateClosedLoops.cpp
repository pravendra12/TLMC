#include "GenerateClosedLoops.h"

vector<vector<size_t>> GenerateClosedLoops(
    const vector<vector<size_t>> &nbrs,
    size_t maxLen,
    size_t numLoops,
    size_t minLen,
    size_t maxTrials,
    size_t maxWalkSteps)
{
  mt19937 rng(random_device{}());

  size_t N = (size_t)nbrs.size();
  uniform_int_distribution<size_t> unisize_t(0, N - 1);

  vector<vector<size_t>> loops;
  loops.reserve(numLoops);

  for (size_t n = 0; n < numLoops; ++n)
  {
    bool found = false;
    for (size_t t = 0; t < maxTrials && !found; ++t)
    {
      size_t start = unisize_t(rng);
      vector<size_t> path;
      path.push_back(start);
      unordered_map<size_t, size_t> indexOf;
      indexOf[start] = 0;

      for (size_t step = 0; step < maxWalkSteps; ++step)
      {
        const auto &neigh = nbrs[path.back()];
        if (neigh.empty())
          break;

        uniform_int_distribution<size_t> uni(0, (size_t)neigh.size() - 1);
        size_t next = neigh[uni(rng)];

        auto it = indexOf.find(next);
        if (it != indexOf.end())
        {
          // Found a cycle
          size_t startIdx = it->second;
          vector<size_t> cycle(path.begin() + startIdx, path.end());
          cycle.push_back(next); // close
          size_t len = (size_t)cycle.size() - 1;
          if (len >= minLen && len <= maxLen)
          {
            loops.push_back(move(cycle));
            found = true;
            break;
          }
        }
        path.push_back(next);
        indexOf[next] = (size_t)path.size() - 1;
      }
    }
  }

  return loops;
}


void WriteLoopStepConfigs(
    const Config &config,
    const vector<size_t> &loop,
    const string &baseDir,
    size_t loopIndex)
{
  // Create a folder for this loop
  fs::path loopDir = fs::path(baseDir) / ("loop_" + to_string(loopIndex));
  if (!fs::exists(loopDir))
    fs::create_directories(loopDir);

  // Copy the original config
  Config configLocal = config;

  for (size_t step = 0; step < loop.size(); ++step)
  {
    size_t latticeId = loop[step];
    // Mark this lattice site
    configLocal.SetElementOfLattice(latticeId, Element("X"));

    // Construct filename
    string filename = (loopDir / ("length_" + to_string(loop.size()) + "step_" + to_string(step) + ".config.gz")).string();

    // Write configuration for this step
    Config::WriteConfig(filename, configLocal);
  }
}