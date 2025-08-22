/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

#include "TimeTemperatureInterpolator.h"

static vector<pair<double, double>> ReadTimeTemperatureJsonFile(const string &filename)
{
  ifstream file(filename);
  if (!file.is_open())
  {
    throw runtime_error("Cannot open file " + filename);
  }

  json jsonData;
  file >> jsonData; // Parse the JSON content into a JSON object

  vector<pair<double, double>> time_temp_vector;
  for (const auto &pair : jsonData)
  {
    double time = pair[0];
    double temp = pair[1];

    time_temp_vector.emplace_back(time, temp);
  }

  return time_temp_vector;
}

TimeTemperatureInterpolator::TimeTemperatureInterpolator(
    const string &timeTemperatureFilename)
{
  if (timeTemperatureFilename.empty())
  {
    return;
  }
  points_ = ReadTimeTemperatureJsonFile(timeTemperatureFilename);
  SortPoints();
}

TimeTemperatureInterpolator::TimeTemperatureInterpolator(
    const vector<pair<double, double>> &points)
    : points_(points)
{
  SortPoints();
}

void TimeTemperatureInterpolator::PrintTimeTemperatureProfile()
{
  for (auto &point : points_)
  {
    cout << point.first << '\t' << point.second << endl;
  }
}

void TimeTemperatureInterpolator::SortPoints()
{
  // Defensive programming. Assume the caller has not sorted the table ascending order
  sort(points_.begin(), points_.end());

  // Ensure that no 2 adjacent x values are equal, atlest we try to divide by zero when we interpolate.
  constexpr double EPSILON{1.0e-12};
  for (size_t i = 1; i < points_.size(); ++i)
  {
    const double deltaX{abs(points_[i].first - points_[i - 1].first)};
    if (deltaX < EPSILON)
    {
      throw range_error("2 adjacent x values are equal." + to_string(points_[i].first));
    }
  }
}

double TimeTemperatureInterpolator::GetTemperature(const double time) const
{
  // Define a lambda that returns true if the time value of a point pair is < the caller's time value
  auto less_than = [](const pair<double, double> &point, double x)
  {
    return point.first < x;
  };

  // Find the first table entry whose value is >= caller's time value
  const auto iter = lower_bound(points_.cbegin(), points_.cend(), time, less_than);

  // If the caller's X value is greater than the largest X value in the table, we can't interpolate.
  if (iter == points_.cend())
  {
    return (points_.cend() - 1)->second;
  }

  // If the caller's X value is less than the smallest X value in the table, we can't interpolate.
  if (iter == points_.cbegin() and time <= points_.cbegin()->first)
  {
    return points_.cbegin()->second;
  }

  // We can interpolate!
  // Equation of line in 2 variables

  const double upper_x{iter->first};
  const double upper_y{iter->second};
  const double lower_x{(iter - 1)->first};
  const double lower_y{(iter - 1)->second};

  const double deltaY{upper_y - lower_y};
  const double deltaX{upper_x - lower_x};
  return lower_y + ((time - lower_x) / deltaX) * deltaY;
}
