#ifndef LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_
#define LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_

#include <omp.h>
#include <filesystem>
#include <nlohmann/json.hpp>
#include "ShortRangeOrder.h"
#include "B2OrderParameter.h"
#include "ClusterDynamics.h"
#include "LocalEnvironment.h"

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
             bool extract_local_env,
             const std::vector<double> &cutoffs_LCE,
             const size_t max_bond_order_LCE,
             const size_t max_cluster_size_LCE);
    virtual ~Traverse();
    void RunAnsys() const;
    //  void RunReformat() const;

  private:
    void RunAnsysOnConfig(
        const size_t configId,
        const Config &config,
        const set<Element> &element_set,
        ostringstream &oss) const;

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

    // Local environment
    const bool extract_LCE_;
    const std::vector<double> cutoffs_LCE_;
    const size_t max_bond_order_LCE_;
    const size_t max_cluster_size_LCE_;
  };

} // namespace ansys

#endif // LMC_LMC_ANSYS_INCLUDE_TRAVERSE_H_