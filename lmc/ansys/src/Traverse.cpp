#include "Traverse.h"

namespace ansys
{
  Traverse::Traverse(
      unsigned long long int initialSteps,
      unsigned long long int incrementSteps,
      const string logType) : initialSteps_(initialSteps),
                              incrementSteps_(incrementSteps),
                              finalSteps_(incrementSteps),
                              logType_(logType),
                              logMap_{},
                              frameOfs_("ansys_log.txt", ofstream::out)
  {

    string logFileName;
    if (logType_ == "kinetic_mc")
    {
      logFileName = "kmc_log.txt";
    }
    else if (logType_ == "canonical_mc")
    {
      logFileName = "cmc_log.txt";
    }
    else if (logType_ == "simulated_annealing")
    {
      logFileName = "sa_log.txt";
    }
    else
    {
      throw(invalid_argument("Unknown log type: " + logType_));
    }

    ifstream ifs(logFileName, ifstream::in);
    if (!ifs.is_open())
    {
      throw(runtime_error("Cannot open " + logFileName));
    }
    while (ifs.peek() == '#')
    {
      ifs.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    string buffer;
    // read header
    getline(ifs, buffer);
    vector<string> headers;
    boost::algorithm::split(headers, buffer, boost::is_any_of("\t"));

    // read data
    while (getline(ifs, buffer))
    {
      if (buffer.empty())
      {
        continue;
      }
      if (buffer[0] == '#')
      {
        continue;
      }
      istringstream line_stream(buffer);
      unsigned long long stepNumber;
      line_stream >> stepNumber;
      if (stepNumber < initialSteps_ || (stepNumber - initialSteps_) % incrementSteps_ != 0)
      {
        continue;
      }
      finalSteps_ = stepNumber;

      unordered_set<string> columns_to_read = {"time", "temperature", "energy"};

      size_t col_index = 1;
      while (line_stream >> buffer)
      {
        const auto &key = headers[col_index];

        if (columns_to_read.find(key) != columns_to_read.end())
        {
          try
          {
            const auto double_value = boost::lexical_cast<double>(buffer);
            if (!holds_alternative<unordered_map<unsigned long long, double>>(logMap_[key]))
            {
              logMap_[key] = unordered_map<unsigned long long, double>();
            }
            get<unordered_map<unsigned long long, double>>(logMap_[key])[stepNumber] = double_value;
          }
          catch (const boost::bad_lexical_cast &)
          {
            if (!holds_alternative<unordered_map<unsigned long long, string>>(logMap_[key]))
            {
              logMap_[key] = unordered_map<unsigned long long, string>();
            }
            get<unordered_map<unsigned long long, string>>(logMap_[key])[stepNumber] = buffer;
          }
        }
        col_index++;
      }
    }
    cout << "Initial Steps: " << initialSteps_ << endl;
    cout << "Increment Steps: " << incrementSteps_ << endl;
    cout << "Final Steps: " << finalSteps_ << endl;

#pragma omp parallel default(none) shared(cout)
    {
#pragma omp master
      {
        cout << "Using " << omp_get_num_threads() << " threads." << endl;
      }
    }
  }

  Traverse::~Traverse() = default;

  void Traverse::RunAnsys(
      const TiledSupercell &tiledSupercell,
      const SubLatticeOccupancy &subLatticeOccupancy,
      const set<Element> &elementSet,
      const unordered_set<size_t> &convertToConfigSet) const
  {
    string outputPath = "";
    if (!convertToConfigSet.empty())
    {
      fs::path cwd = fs::current_path();
      fs::path configDir = cwd / "config";

      if (!fs::exists(configDir))
        fs::create_directory(configDir);

      outputPath = configDir.string();

      std::cout << "Config directory created at: " << outputPath << std::endl;
    }

    // Write header once
    frameOfs_ << GetHeaderFrameString(elementSet) << flush;

    // Generate config indices
    for (unsigned long long i = initialSteps_; i <= finalSteps_; i += incrementSteps_)
    {
      // Read atomic indices (or atomic numbers) from the compressed binary file
      const string atomicIndicesFilename = to_string(i) + ".bin.gz";

      if (!fs::exists(atomicIndicesFilename))
      {
        continue;
      }

      cout << "Running Ansys for " << atomicIndicesFilename << endl;

      auto atomicIndicesVector = TiledSupercell::ReadAtomicIndicesFromFile(atomicIndicesFilename);

      const auto time = logMap_.find("time") == logMap_.end()
                            ? nan("")
                            : get<unordered_map<unsigned long long, double>>(logMap_.at("time")).at(i);

      const auto temperature = get<unordered_map<unsigned long long, double>>(logMap_.at("temperature")).at(i);
      const auto energy = get<unordered_map<unsigned long long, double>>(logMap_.at("energy")).at(i);

      frameOfs_ << i << "\t" << time << "\t" << temperature << "\t" << energy;

      bool saveConfig = false;
      if (convertToConfigSet.find(i) != convertToConfigSet.end())
      {
        saveConfig = true;
      }

      string filename = "";

      if (!outputPath.empty())
      {
        filename = outputPath + "/" + to_string(i) + ".xyz.gz";
      }

      ostringstream oss;

      RunAnsysOnConfig(
          tiledSupercell,
          atomicIndicesVector,
          subLatticeOccupancy,
          elementSet,
          oss,
          saveConfig,
          filename);

      frameOfs_ << oss.str() << "\n";
    }

    frameOfs_.close();
  }

  void Traverse::RunAnsysOnConfig(
      const TiledSupercell &tiledSupercell,
      const vector<size_t> &atomicIndicesVector,
      const SubLatticeOccupancy &subLatticeOccupancy,
      const set<Element> &elementSet,
      ostringstream &oss,
      const bool &saveConfig,
      const string &outfilename) const
  {
    // Analysis

    /// B2 Order Parameter
    const auto elementOccupanciesMap = subLatticeOccupancy.GetAlphaBetaSiteOccupancy(
        atomicIndicesVector);

    const auto orderParamMap = subLatticeOccupancy.ComputeOrderParameter(
        elementOccupanciesMap);

    for (const auto &elementOccupancy : elementOccupanciesMap)
    {
      auto element = elementOccupancy.first;

      auto alphaOccupancy = elementOccupanciesMap.at(element).first;
      auto betaOccupancy = elementOccupanciesMap.at(element).second;
      auto orderParam = orderParamMap.at(element);

      if (elementSet.find(element) != elementSet.end())
      {
        oss << "\t" << orderParam << "\t" << alphaOccupancy << "\t" << betaOccupancy;
      }
    }

    /// WCP Parameter
    auto tiledSupercellLocal = tiledSupercell;
    tiledSupercellLocal.UpdateAtomVector(atomicIndicesVector);

    ShortRangeOrderTLMC sroObject(tiledSupercellLocal);

    auto sroParamMap1 = sroObject.ComputeWarrenCowley(1);
    auto sroParamMap2 = sroObject.ComputeWarrenCowley(2);

    for (const auto &pair : sroParamMap1)
    {
      double sro1Value = sroParamMap1.count(pair.first) ? sroParamMap1.at(pair.first) : nan("");
      double sro2Value = sroParamMap2.count(pair.first) ? sroParamMap2.at(pair.first) : nan("");
      oss << "\t" << sro1Value << "\t" << sro2Value;
    }

    if (saveConfig)
    {
      B2ClusterTLMC b2Cluster(tiledSupercellLocal);
      b2Cluster.WriteB2ClusterConfig(outfilename);
    }
  }

  string Traverse::GetHeaderFrameString(const set<Element> &elementSet) const
  {
    string headerFrame_ = "steps\ttime\ttemperature\tenergy\t";

    // B2 Order Parameter
    for (const auto &element : elementSet)
    {
      auto elementString = element.GetElementString();
      headerFrame_ += "B2_order_param_" + elementString + "\t" +
                      "alpha_occupancy_" + elementString + "\t" +
                      "beta_occupancy_" + elementString + "\t";
    }

    // WCP Short range order parameter
    // Right now the order list is hard coded as we would update the neighbourlist
    // upto second nn for getting the B2Clusters information

    static const vector<string> order_list{"first", "second"};
    for (auto element1 : elementSet)
    {
      for (auto element2 : elementSet)
      {
        for (const auto &order : order_list)
        {
          headerFrame_ += "warren_cowley_" + order + "_" + element1.GetElementString() + "-" + element2.GetElementString() + "\t";
        }
      }
    }

    if (!headerFrame_.empty() && headerFrame_.back() == '\t')
    {
      headerFrame_.back() = '\n';
    }

    return headerFrame_;
  }

} // namespace ansys