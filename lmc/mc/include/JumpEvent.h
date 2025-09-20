/*******************************************************************************
 * Copyright (c) 2022-2025. All rights reserved.
 * @Author: Zhucong Xi
 * @Date: 2022
 * @Last Modified by: pravendra12
 * @Last Modified: 2025-06-01
 ******************************************************************************/

/*! \file JumpEvent.h
 *  @brief File for JumpEvent class declaration.
 */

#ifndef LMC_MC_INCLUDE_JUMPEVENT_H_
#define LMC_MC_INCLUDE_JUMPEVENT_H_

#include <cstddef>
#include <utility>
#include <array>
#include <cmath>
#include "LatticeSiteMapping.hpp"

using namespace std;

namespace mc
{
  class JumpEvent
  {
  public:
    /** @brief Default Constructor
     */
    JumpEvent();

    /** @brief Value Constructor
     *  @param jumpPair Lattice Jump Pair
     *  @param barrierAndDe Barrier and Energy change
     *  @param beta Inverse temperature factor (1/kT), used in rate calculations.
     */
    JumpEvent(
        const pair<LatticeSiteMapping, LatticeSiteMapping> &jumpPair,
        const pair<double, double> &barrierAndDe,
        double beta);

    /**
     * @brief Get the Id Jump Pair
     */
    [[nodiscard]] const pair<LatticeSiteMapping, LatticeSiteMapping> &GetIdJumpPair() const;

    /**
     * @brief Get the Forward Barrier
     */
    [[nodiscard]] double GetForwardBarrier() const;

    /**
     * @brief Get the Forward Rate
     */
    [[nodiscard]] double GetForwardRate() const;

    /**
     * @brief Get the Backward Barrier
     */
    [[nodiscard]] double GetBackwardBarrier() const;

    /**
     * @brief Get the Backward Rate
     */
    [[nodiscard]] double GetBackwardRate() const;

    /**
     * @brief Get the Energy Change
     */
    [[nodiscard]] double GetEnergyChange() const;

    /**
     * @brief Get the Probability
     */
    [[nodiscard]] double GetProbability() const;

    /**
     * @brief Get the Cumulative Probability
     */
    [[nodiscard]] double GetCumulativeProbability() const;

    /**
     * @brief Get the Reverse Jump Event
     */
    [[nodiscard]] JumpEvent GetReverseJumpEvent() const;

    /**
     * @brief Set the Probability
     *
     * @param probability Probability of jump
     */
    void SetProbability(double probability);

    /**
     * @brief Set the Cumulative Probability object
     *
     * @param cumulative_probability Cumulative Probability
     */
    void SetCumulativeProbability(double cumulative_probability);

    /**
     * @brief Calculate probabiltiy of the jump
     *
     * @param total_rates Total rate
     */
    void CalculateProbability(double total_rates);

  private:
    /**
     * @brief Inverse temperature factor (1/kT), used in rate calculations.
     */
    double beta_{};

    /**
     * @brief Pair of lattice site indices representing the atom-vacancy jump.
     */
    std::pair<LatticeSiteMapping, LatticeSiteMapping> jump_pair_{};

    /**
     * @brief Energy barrier for the forward jump (in eV).
     */
    double forward_barrier_{};

    /**
     * @brief Energy change due to the jump event (in eV).
     */
    double energy_change_{};

    /**
     * @brief Rate of the forward jump computed using the Arrhenius expression.
     */
    double forward_rate_{};

    /**
     * @brief Probability of selecting this jump event among all possible events.
     */
    double probability_{};

    /**
     * @brief Cumulative probability used in event selection (e.g., in rejection-free kMC).
     */
    double cumulative_probability_{};
  };
} // mc
#endif // LMC_MC_INCLUDE_JUMPEVENT_H_
