#include "Config.h"
#include <vector>
#include <unordered_map>
#include <random>
#include <iostream>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

vector<vector<size_t>> GenerateClosedLoops(
    const vector<vector<size_t>> &nbrs,
    size_t maxLen,
    size_t numLoops,
    size_t minLen = 3,
    size_t maxTrials = 5000,
    size_t maxWalkSteps = 2000);

void WriteLoopStepConfigs(
    const Config &config,
    const vector<size_t> &loop,
    const string &baseDir,
    size_t loopIndex);
