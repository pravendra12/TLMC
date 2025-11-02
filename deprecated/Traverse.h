#ifndef LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_
#define LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_

#include <omp.h>
#include <filesystem>
#include <nlohmann/json.hpp>
// #include "ShortRangeOrder.h"
#include "B2OrderParameter.h"
// #include "ClusterDynamics.h"
// #include "ConfigEncoding.h"
#include <memory>

namespace fs = std::filesystem;

namespace ansys
{

  class Traverse
  {
  public:
    Traverse(unsigned long long int initial_steps,
             unsigned long long int increment_steps,
             const std::vector<double> &cutoffs,
             std::string log_type,
             std::string config_type,
             const bool extract_encoding,
             const size_t maxBondOrder, 
             const size_t maxClusterSize);
    virtual ~Traverse();
    void RunAnsys() const;
    //  void RunReformat() const;


    static void RunAnsysOnConfig(
      const size_t configId, 
      const Config &config, 
      const set<Element> &element_set, 
      ostringstream &oss, 
      const string &outputFolder);
    
    private:

    void RunAnsysLCE(
        const size_t configId,
        const Config &config,
        set<Element> element_set,
        ostringstream &oss) const;

    std::string GetHeaderFrameString(const std::set<Element> &element_set) const;

    std::string GetHeaderFrameStringWithFlag(const std::set<Element> &element_set, const string &flag) const;

    const unsigned long long initial_steps_;
    const unsigned long long increment_steps_;
    unsigned long long final_steps_;
    std::vector<double> cutoffs_;
    const std::string log_type_;
    const std::string config_type_;

    using MapVariant =
        std::variant<std::unordered_map<unsigned long long, double>, std::unordered_map<unsigned long long, std::string>>;

    std::unordered_map<std::string, MapVariant> log_map_;

    mutable std::ofstream frame_ofs_;
    
    // Encoding
    const bool extract_encoding_{};
    const size_t maxBondOrder_{};
    const size_t maxClusterSize_{};

    void RunAnsys(const B2OrderParameter &b2OrderParam, const set<Element> &elementSet) const;
  
  };

} // namespace ansys

#endif // LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_