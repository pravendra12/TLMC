/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/22/23 10:40 PM                                                          *
 **************************************************************************************************/

#include "ShortRangeOrder.h"

#include <utility>

ShortRangeOrder::ShortRangeOrder(
    const Config &config,
    const std::set<Element> &element_set) : config_(config),
                                            element_set_(element_set),
                                            concentration_(config_.GetConcentration())
{
}

std::map<std::string, double> ShortRangeOrder::FindWarrenCowley(const size_t shell_number) const
{
  std::map<std::string, double> warren_cowley;

  // Initializing the warren_cowley map
  for (const auto &element1 : element_set_)
  {
    for (const auto &element2 : element_set_)
    {
      if (element1 == ElementName::X || element2 == ElementName::X)
      {
        continue;
      }
      warren_cowley[element1.GetElementString() + "-" + element2.GetElementString()] = 0;
    }
  }

  auto element_list_map = config_.GetElementOfAtomIdVectorMap();

  auto num_bonds = config_.GetNeighborLatticeIdVectorOfLattice(0, shell_number).size();

  for (const auto &element1 : element_set_)
  {

    if (element1 == ElementName::X)
    {
      continue;
    }
    // size_t num_all_bonds = element_list_map.at(element1).size() * num_bonds;

    size_t num_all_bonds = element_list_map.at(element1).size() * num_bonds;

    std::map<Element, size_t> ct_this_pair_map{};
    for (const auto &element2 : element_set_)
    {
      if (element2 == ElementName::X)
      {
        continue;
      }
      ct_this_pair_map[element2] = 0;
    }
    for (const auto &atom_id1 : element_list_map.at(element1))
    {
      std::vector<size_t> neighbor_list = config_.GetNeighborAtomIdVectorOfAtom(atom_id1, shell_number);

      for (const auto &atom_id2 : neighbor_list)
      {
        const auto &this_element = config_.GetElementOfAtom(atom_id2);
        if (this_element == ElementName::X)
        {
          continue;
        }
        ct_this_pair_map[this_element]++;
      }
    }
    for (auto [element2, ct_this_pair] : ct_this_pair_map)
    {
      std::string key = element1.GetElementString() + "-" + element2.GetElementString();
      double pij = static_cast<double>(ct_this_pair) / static_cast<double>(num_all_bonds);
      warren_cowley[key] = (pij - concentration_.at(element2)) /
                           (static_cast<double>(element1 == element2) - concentration_.at(element2));
      // res[key] = 1 - (pij / concentration_[element2]);
    }
  }
  return warren_cowley;
}
