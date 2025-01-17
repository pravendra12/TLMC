/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: pravendra12                                                                 *
 * @Last Modified time: 11/30/24 5:50 PM                                                          *
 **************************************************************************************************/

#ifndef LMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_
#define LMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_
#include <list>
#include <deque>
#include <iterator>

/*! \file ThermodynamicAveraging.h
 *  \brief  File for ThermodynamicAveraging class declaration.
 */
namespace mc {

class ThermodynamicAveraging {
 public:

  /*! \brief Constructs a ThermodynamicAveraging object.
   *  \param size The maximum size of the sliding window for averaging.
   */
  explicit ThermodynamicAveraging(size_t size);

  /*! \brief Adds a new energy value.
   *  \param value The energy value to be added.
   */
  void AddEnergy(double value);

  /*! \brief Computes the thermodynamic average for the current sliding window.
   * 
   *  This method calculates the weighted average of energy values using a 
   *  Boltzmann-like distribution. It incorporates the inverse temperature 
   *  (`beta`) and normalizes by the partition function.
   * 
   *  \param beta The inverse temperature (1/kT) for weighting the average.
   *  \return The computed thermodynamic average.
   */
  [[nodiscard]] double GetThermodynamicAverage(double beta) const;

 private:

  /*! \brief The arithmetic average of the energy values in the sliding window.
   * 
   *  If the window is empty, the average is returned as 0.
   * 
   *  \return The arithmetic average of the energy values.
   */
  [[nodiscard]] double GetAverage() const;

  /// Energy List.
  std::deque<double> energy_list_{};

  /// The maximum size of the sliding window.
  const size_t size_;

  /// The sum of the energy values in the sliding window.
  double sum_{};
};

} // mc

#endif //LMC_MC_INCLUDE_THERMODYNAMICAVERAGING_H_
