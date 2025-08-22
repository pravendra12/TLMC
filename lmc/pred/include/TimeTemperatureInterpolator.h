/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

#ifndef LMC_PRED_INCLUDE_TIMETEMPERATUREINTERPOLATOR_H_
#define LMC_PRED_INCLUDE_TIMETEMPERATUREINTERPOLATOR_H_

/*! @file TimeTemperatureInterpolator.h
    @brief Interpolates temperature based on the cooling curve
 */

#include <iostream>
#include <utility>
#include <vector>
#include <algorithm>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
using namespace std;

class TimeTemperatureInterpolator
{
public:
  /*! @brief Constructor for time temperature interpolator
      @param timeTemperatureFilename Time temperature json file
  */
  explicit TimeTemperatureInterpolator(const string &timeTemperatureFilename);

  /*! @brief Constructor for time temperature interpolator
      @param points (Time-Temperature) Vector
  */
  explicit TimeTemperatureInterpolator(const vector<pair<double, double>> &points);

  /*! @brief Computes the corresponding Y value for X using linear interpolation
   */
  [[nodiscard]] double GetTemperature(double time) const;

  /*! @brief Prints the time temperature profile
   */
  void PrintTimeTemperatureProfile();

private:
  /*! @brief Sort the points in ascending order based on time
   */
  void SortPoints();

  /*! @brief Contains Time-Temperature pairs
   */
  vector<pair<double, double>> points_{};
};

#endif // LMC_PRED_INCLUDE_TIMETEMPERATUREINTERPOLATOR_H_
