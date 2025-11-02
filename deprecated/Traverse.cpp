#include "Traverse.h"

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/lexical_cast.hpp>
#include <iostream>
#include <omp.h>

namespace ansys
{
  static Config GetConfig(
      const std::string &config_type, size_t i, std::vector<double> cutoff)
  {
    Config config;
    if (config_type == "config")
    {
      string base = std::to_string(i);
      try
      {
        config = Config::ReadConfig(base + ".cfg");
      }
      catch (...)
      {
        try
        {
          config = Config::ReadConfig(base + ".cfg.gz");
        }
        catch (...)
        {
          throw std::runtime_error("Failed to load config: tried " + base + ".cfg and " + base + ".cfg.gz");
        }
      }
      config.UpdateNeighborList(cutoff);
    }
    else
    {
      throw std::invalid_argument("Unknown config type: " + config_type);
    }
    return config;
  }

  Traverse::Traverse(unsigned long long int initial_steps,
                     unsigned long long int increment_steps,
                     const std::vector<double> &cutoffs,
                     std::string log_type,
                     std::string config_type,
                     const bool extract_encoding,
                     const size_t maxBondOrder,
                     const size_t maxClusterSize)
      : initial_steps_(initial_steps),
        increment_steps_(increment_steps),
        final_steps_(increment_steps),
        cutoffs_(std::move(cutoffs)),
        log_type_(std::move(log_type)),
        config_type_(std::move(config_type)),
        log_map_{},
        frame_ofs_("ansys_frame_log.txt", std::ofstream::out),
        extract_encoding_(extract_encoding),
        maxBondOrder_(maxBondOrder),
        maxClusterSize_(maxClusterSize)
  {

    std::string log_file_name;
    if (log_type_ == "kinetic_mc")
    {
      log_file_name = "kmc_log.txt";
    }
    else if (log_type_ == "canonical_mc")
    {
      log_file_name = "cmc_log.txt";
    }
    else if (log_type_ == "simulated_annealing")
    {
      log_file_name = "sa_log.txt";
    }
    else
    {
      throw std::invalid_argument("Unknown log type: " + log_type_);
    }

    std::ifstream ifs(log_file_name, std::ifstream::in);
    if (!ifs.is_open())
    {
      throw std::runtime_error("Cannot open " + log_file_name);
    }
    while (ifs.peek() == '#')
    {
      ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }

    std::string buffer;
    // read header
    std::getline(ifs, buffer);
    std::vector<std::string> headers;
    boost::algorithm::split(headers, buffer, boost::is_any_of("\t"));

    // read data
    while (std::getline(ifs, buffer))
    {
      if (buffer.empty())
      {
        continue;
      }
      if (buffer[0] == '#')
      {
        continue;
      }
      std::istringstream line_stream(buffer);
      unsigned long long step_number;
      line_stream >> step_number;
      if (step_number < initial_steps_ || (step_number - initial_steps_) % increment_steps != 0)
      {
        continue;
      }
      final_steps_ = step_number;
      size_t col_index = 1;
      while (line_stream >> buffer)
      {
        const auto &key = headers[col_index];
        try
        {
          const auto double_value = boost::lexical_cast<double>(buffer);
          if (!std::holds_alternative<std::unordered_map<unsigned long long, double>>(log_map_[key]))
          {
            log_map_[key] = std::unordered_map<unsigned long long, double>();
          }
          std::get<std::unordered_map<unsigned long long, double>>(log_map_[key])[step_number] = double_value;
        }
        catch (const boost::bad_lexical_cast &)
        {
          if (!std::holds_alternative<std::unordered_map<unsigned long long, std::string>>(log_map_[key]))
          {
            log_map_[key] = std::unordered_map<unsigned long long, std::string>();
          }
          std::get<std::unordered_map<unsigned long long, std::string>>(log_map_[key])[step_number] = buffer;
        }
        col_index++;
      }
    }

    std::cout << "Initial Steps: " << initial_steps_ << std::endl;
    std::cout << "Increment Steps: " << increment_steps_ << std::endl;
    std::cout << "Final Steps: " << final_steps_ << std::endl;

#pragma omp parallel default(none) shared(std::cout)
    {
#pragma omp master
      {
        std::cout << "Using " << omp_get_num_threads() << " threads." << std::endl;
      }
    }
  }

  Traverse::~Traverse() = default;

  // Parallel Version
  void Traverse::RunAnsys() const
  {
    std::set<Element> element_set;
    std::vector<unsigned long long> config_indices;

    // Generate config indices
    for (unsigned long long i = initial_steps_; i <= final_steps_; i += increment_steps_)
    {
      config_indices.push_back(i);
    }

    int total_configs = static_cast<int>(config_indices.size());
    int num_threads = omp_get_max_threads();

    const Config referenceConfig;
    const set<Element> elementSet;

    std::unique_ptr<ConfigEncoding> configEncoder;

    // Write header once
    {
      auto config = GetConfig(config_type_, config_indices[0], cutoffs_);
      auto atomVector = config.GetAtomVector();
      element_set = std::set<Element>(atomVector.begin(), atomVector.end());

      if (extract_encoding_)
      {
        // Declare configEncoder once
        configEncoder = std::make_unique<ConfigEncoding>(
            config,
            element_set,
            maxBondOrder_,
            maxClusterSize_);
      }

      element_set.erase(Element("X"));
      frame_ofs_ << GetHeaderFrameString(element_set) << std::flush;
    }

    // Process in batches
    for (int batch_start = 0; batch_start < total_configs; batch_start += num_threads)
    {
      int batch_end = std::min(batch_start + num_threads, total_configs);
      std::vector<std::ostringstream> output_buffers(batch_end - batch_start);

// Parallel region
#pragma omp parallel for num_threads(num_threads)
      for (int idx = batch_start; idx < batch_end; ++idx)
      {
        int local_index = idx - batch_start;
        unsigned long long i = config_indices[idx];
        auto config = GetConfig(config_type_, i, cutoffs_);

        const auto time = log_map_.find("time") == log_map_.end()
                              ? nan("")
                              : std::get<std::unordered_map<unsigned long long, double>>(log_map_.at("time")).at(i);

        const auto average_time = log_map_.find("average_time") == log_map_.end()
                                      ? nan("")
                                      : std::get<std::unordered_map<unsigned long long, double>>(log_map_.at("average_time")).at(i);

        const auto temperature = std::get<std::unordered_map<unsigned long long, double>>(log_map_.at("temperature")).at(i);
        const auto energy = std::get<std::unordered_map<unsigned long long, double>>(log_map_.at("energy")).at(i);

        size_t selected_lattice_id = SIZE_MAX;

        if (log_type_ == "kinetic_mc")
        {
          selected_lattice_id = static_cast<size_t>(std::get<std::unordered_map<unsigned long long, double>>(log_map_.at("selected")).at(i));
        }

        std::ostringstream &oss = output_buffers[local_index];
        oss << i << "\t" << time << "\t" << average_time << "\t" << temperature << "\t" << energy;

        // Analysis on the original configuration
        RunAnsysOnConfig(i, config, element_set, oss, "AnalyzedConfigs");

        if (extract_encoding_)
        {
          // Write the ce encoding to the file
          oss << "\t";
          VectorXd encodingVectorBefore = configEncoder->GetEncodeVector(config);

          oss << encodingVectorBefore.transpose();

          // Lattice Jump
          oss << "\t";
          size_t vacancyId = config.GetVacancyLatticeId();

          // cout << "VacancyId: " << vacancyId << endl;
          // cout << "Selected_Id: " << selected_lattice_id << endl;
          // cout << "Distance Order: " << config.GetDistanceOrder(vacancyId, selected_lattice_id) << endl;

          config.LatticeJump(make_pair(vacancyId, selected_lattice_id));
          VectorXd encodingVectorAfter = configEncoder->GetEncodeVector(config);


          // Also need to run the analysis on this config 
          // but in principle these two config are just one step away hence we can assume
          // that the sro param and B2 order param will be nearly same

          oss << encodingVectorAfter.transpose();
        }

        oss << "\n";
      }

      // Write the batch to file in order
      for (auto &oss : output_buffers)
      {
        frame_ofs_ << oss.str();
      }
    }

    frame_ofs_.close();
  }

  void Traverse::RunAnsysOnConfig(
      const size_t configId,
      const Config &config,
      const set<Element> &element_set,
      ostringstream &oss,
      const string &outputFolder)
  {
    // Analysis

    /// Short range order

    ShortRangeOrder short_range_order(config, element_set);
    const auto sro1 = short_range_order.FindWarrenCowley(1);
    const auto sro2 = short_range_order.FindWarrenCowley(2);
    const auto sro3 = short_range_order.FindWarrenCowley(3);

    for (const auto &pair : sro1)
    {
      double sro1_value = sro1.count(pair.first) ? sro1.at(pair.first) : nan("");
      double sro2_value = sro2.count(pair.first) ? sro2.at(pair.first) : nan("");
      double sro3_value = sro3.count(pair.first) ? sro3.at(pair.first) : nan("");
      oss << "\t" << sro1_value << "\t" << sro2_value << "\t" << sro3_value;
    }

    // Global List
    map<string, Config::ValueVariant> globalList;

    /// B2 Order

    B2OrderParameter b2Order(config);
    for (auto element : element_set)
    {
      double b2OrderParameter = b2Order.GetOrderParameter(element);
      double alphaOccupancy = b2Order.GetAlphaSiteOccupancy(element);
      double betaOccupancy = b2Order.GetBetaSiteOccupancy(element);

      globalList["b2OrderParameter" + element.GetElementString()] = b2OrderParameter;

      oss << "\t" << b2OrderParameter << "\t" << alphaOccupancy << "\t" << betaOccupancy;
    }

    // Cluster Dynamics

    // Auxilary List
    std::map<std::string, Config::VectorVariant> auxiliaryList;

    // Cluster Size
    vector<int> clusterSizeVector;

    ClusterDynamics b2Cluster(config);

    b2Cluster.detectB2Clusters(auxiliaryList, clusterSizeVector);

    // Write cluster size to log file
    oss << "\t";
    for (size_t i = 0; i < clusterSizeVector.size(); ++i)
    {
      if (i == clusterSizeVector.size() - 1)
      {
        oss << clusterSizeVector[i];
      }
      else
      {
        oss << clusterSizeVector[i] << ", ";
      }
    }

    if (!fs::exists(outputFolder))
    {
      fs::create_directories(outputFolder);
    }

    // Write to file
    Config::WriteXyzExtended(outputFolder + "/" + to_string(configId) + ".xyz.gz",
                             config,
                             auxiliaryList,
                             globalList);
  }

  std::string Traverse::GetHeaderFrameString(const std::set<Element> &element_set) const
  {
    std::string header_frame = "steps\ttime\taverage_time\ttemperature\tenergy\t";

    // For simple analysis
    header_frame += GetHeaderFrameStringWithFlag(element_set, "");

    // Local Config
    if (extract_encoding_)
    {
      header_frame += "ce_encoding_before_jump\tce_encoding_after_jump\t";
    }

    if (!header_frame.empty() && header_frame.back() == '\t')
    {
      header_frame.back() = '\n';
    }

    return header_frame;
  }

  std::string Traverse::GetHeaderFrameStringWithFlag(const std::set<Element> &element_set, const string &flag) const
  {
    std::string header_frame = "";
    // SRO Parameter
    static const std::vector<std::string> order_list{"first", "second", "third"};
    for (auto element1 : element_set)
    {
      for (auto element2 : element_set)
      {
        for (const auto &order : order_list)
        {
          header_frame += "warren_cowley_" + order + "_" + element1.GetElementString() + "-" + element2.GetElementString() + flag + "\t";
        }
      }
    }

    // B2 Order Parameter

    for (auto element : element_set)
    {
      auto elementString = element.GetElementString();
      header_frame += "B2_order_param_" + elementString + flag + "\t" +
                      "alpha_occupancy_" + elementString + flag + "\t" +
                      "beta_occupancy_" + elementString + flag + "\t";
    }

    // Cluster Dynamics

    header_frame += "B2_cluster_size" + flag + "\t";

    return header_frame;
  }

} // namespace ansys