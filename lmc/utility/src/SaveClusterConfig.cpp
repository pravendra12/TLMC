#include "SaveClusterConfig.h"

void SaveClustersConfig(const string outputPath,
                        const Config &config,
                        const vector<set<vector<size_t>>> &equivalentClusters)
{

  int i = 0;
  for (const auto &orbit : equivalentClusters)
  {
    cout << "Orbit Index : " << i << endl;

    // Get lattice cluster type (for naming)
    auto latticeClusterType = IdentifyLatticeClusterType(config, *orbit.begin());
    cout << latticeClusterType << endl;

    // First, create a copy of config where ALL clusters in the orbit are set to Oxygen
    auto configOxygen = config;
    for (const auto &cluster : orbit)
    {
      for (auto id : cluster)
      {
        configOxygen.SetElementOfLattice(id, Element("O"));
      }
    }

    // Now, iterate over each cluster in the orbit
    int j = 0;
    for (const auto &cluster : orbit)
    {
      auto configCluster = configOxygen; // start from Oxygen-modified config

      // Set this cluster's sites to Nitrogen
      for (auto id : cluster)
      {
        configCluster.SetElementOfLattice(id, Element("N"));
      }

      // Build filename
      ostringstream oss;
      oss << "Orbit_" << i << "_" << latticeClusterType
          << "_cluster_" << j << ".cfg.gz";
      auto filename = outputPath + "/" + oss.str();

      cout << "Writing: " << filename << endl;

      // Write config
      Config::WriteConfig(filename, configCluster);

      j++;
    }

    i++;
  }
}
