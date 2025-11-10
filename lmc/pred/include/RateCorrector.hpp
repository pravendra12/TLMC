#ifndef LMC_LMC_PRED_INCLUDE_RATECORRECTOR_HPP_
#define LMC_LMC_PRED_INCLUDE_RATECORRECTOR_HPP_
#include <cmath>
#include <utility>
#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Constants.hpp"
#include "Element.hpp"

using namespace std;

class RateCorrector
{

public:
  RateCorrector(
      const map<Element, double> &concentrationMap,
      const string &pathVFEOutput) : concentrationMap_(concentrationMap),
                                     vfeMap_(GetVfeMap(pathVFEOutput))
  {
  }

  [[nodiscard]] inline double GetTimeCorrectionFactor(double temperature)
  {
    auto vacancy = Element("X");
    return concentrationMap_.at(vacancy) / GetCorrectVacancyConcentration(temperature);
  }

private:
  // Read the vfe and the do and ensemble average

  inline double GetCorrectVacancyConcentration(double temperature)
  {
    double avgCv = 0; // average vacancy concentration
    UpdateAverageVacancyConcentration(temperature);

    for (const auto &[element, conc] : concentrationMap_)
    {
      avgCv += conc * vacancyConcentrationMap_[element][temperature];
    }
    return avgCv;
  }

  void UpdateAverageVacancyConcentration(const double temperature)
  {
    const double beta = 1 / constants::kBoltzmann / temperature;

    for (const auto &[eleString, vfeArray] : vfeMap_)
    {
      auto element = Element(eleString);

      if (vacancyConcentrationMap_[element].find(temperature) != vacancyConcentrationMap_[element].end())
        continue;

      double avgConc = 0;
      for (size_t i = 0; i < vfeArray.size(); i++)
        avgConc += std::exp(-beta * vfeArray[i]);

      avgConc /= double(vfeArray.size());

      vacancyConcentrationMap_[element][temperature] = avgConc;
    }
  }

  map<string, vector<double>> GetVfeMap(const string &pathVFEOutput)
  {
    ifstream infile(pathVFEOutput);

    if (!infile.is_open())
    {
      cerr << "Error in `RateCorrector` : Can't open file: " << pathVFEOutput << endl;
    }

    map<string, vector<double>> vfeMap;
    string line;

    // idx  Mo  Ta
    // 0  2.1 3.2
    // ..
    // N  2.3 3.01
    getline(infile, line);
    stringstream ssHeader(line);
    string col;
    vector<string> columns;

    while (getline(ssHeader, col, '\t'))
    {
      if (col != "idx")
      {
        columns.push_back(col);
        vfeMap[col] = {};
      }
    }

    while (getline(infile, line))
    {
      if (line.empty())
        continue;

      stringstream ss(line);
      string item;

      // First item = idx, skip
      getline(ss, item, '\t');

      for (const auto &colName : columns)
      {
        if (!getline(ss, item, '\t'))
          break;
        vfeMap[colName].push_back(stod(item));
      }
    }

    infile.close();

    return vfeMap;
  }

  const map<Element, double> concentrationMap_;
  const map<string, vector<double>> vfeMap_;

  // Element : <Temperature, VacancyConcentration>
  map<Element, map<double, double>> vacancyConcentrationMap_;
};

#endif // LMC_LMC_PRED_INCLUDE_RATECORRECTOR_HPP_
