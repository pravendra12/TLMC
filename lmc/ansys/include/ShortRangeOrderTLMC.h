/**************************************************************************************************
 * Copyright (c) 2022-2023. All rights reserved.                                                  *
 * @Author: Zhucong Xi                                                                            *
 * @Date: 5/30/22 4:05 PM                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 8/22/23 10:40 PM                                                          *
 **************************************************************************************************/

#ifndef LMC_ANSYS_INCLUDE_SHORTRANGEORDERTLMC_H_
#define LMC_ANSYS_INCLUDE_SHORTRANGEORDERTLMC_H_

#include "TiledSupercell.h"
#include <utility>

using namespace std;

using PairCountByCenter = std::unordered_map<std::string, std::unordered_map<std::string, size_t>>;

class ShortRangeOrderTLMC
{
public:
  ShortRangeOrderTLMC(
      const TiledSupercell &tiledSupercell);

  map<string, double> ComputeWarrenCowley(const size_t shellNumber);

  PairCountByCenter GetElementPairCountMap(const size_t shellNumber);

protected:
  PairCountByCenter ComputeElementPairCountsForSmallConfig(
    const size_t smallConfigId, 
    const size_t shellNumber);
    
  string GetElementPairString(
      const Element &element1,
      const Element &element2);

  const TiledSupercell &tiledSupercell_;
  const map<Element, double> compositionMap_;
};

#endif // LMC_ANSYS_INCLUDE_SHORTRANGEORDERTLMC_H_
