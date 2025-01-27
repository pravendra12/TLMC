/**************************************************************************************************
 * Copyright (c) 2023. All rights reserved.                                                       *
 * @Author: Zhucong Xi                                                                            *
 * @Date:                                                                                         *
 * @Last Modified by: zhucongx                                                                    *
 * @Last Modified time: 10/30/23 3:09 PM                                                          *
 **************************************************************************************************/

#include <cmath>
#include <utility>
#include "JumpEvent.h"
namespace mc {
JumpEvent::JumpEvent() = default;
JumpEvent::JumpEvent(std::pair<size_t, size_t> jump_pair,
                     const std::array<double, 3> &barrier_and_diff,
                     double beta)
    : beta_(beta),
      jump_pair_(std::move(jump_pair)),
      forward_barrier_(barrier_and_diff[0]),
      backward_barrier_(barrier_and_diff[1]),
      energy_change_(barrier_and_diff[2]),
      forward_rate_(std::exp(-forward_barrier_ * beta_)) {}

      
const std::pair<size_t, size_t> &JumpEvent::GetIdJumpPair() const {
  return jump_pair_;
}
double JumpEvent::GetForwardBarrier() const {
  return forward_barrier_;
}
double JumpEvent::GetForwardRate() const {
  return forward_rate_;
}
double JumpEvent::GetBackwardBarrier() const {
  return backward_barrier_;
}
double JumpEvent::GetBackwardRate() const {
  return std::exp(-backward_barrier_ * beta_);
}
double JumpEvent::GetEnergyChange() const {
  return energy_change_;
}

double JumpEvent::GetProbability() const {
  return probability_;
}

double JumpEvent::GetCumulativeProbability() const {
  return cumulative_probability_;
}

void JumpEvent::SetProbability(double probability) {
  probability_ = probability;
}

void JumpEvent::SetCumulativeProbability(double cumulative_probability) {
  cumulative_probability_ = cumulative_probability;
}

void JumpEvent::CalculateProbability(double total_rates) {
  probability_ =  forward_rate_/total_rates;
}
JumpEvent JumpEvent::GetReverseJumpEvent() const {
  std::array<double, 3> backward_barrier_and_diff;
  backward_barrier_and_diff[0] = backward_barrier_;
  backward_barrier_and_diff[1] = forward_barrier_;
  backward_barrier_and_diff[2] = -energy_change_;

  return JumpEvent{{jump_pair_.second, jump_pair_.first},
                   backward_barrier_and_diff, beta_};
}
} // mc
