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
#include <numbers>
#include "Constants.hpp"
#include "Element.hpp"

using namespace std;

class RateCorrector
{
  // Current implementation only accounts for the constant temperature vacancy
  // evolution given a initial vacancy concentration along with other parameters.
  // Later temperature dependence can be added, which can be used during quenching
  // simulations.

  /*
  Reference:
    Fischer, F. D., J. Svoboda, F. Appel, and E. Kozeschnik.
    "Modeling of Excess Vacancy Annihilation at Different Types of Sinks."
    Acta Materialia 59(9) (2011) 3463–3472. https://doi.org/10.1016/j.actamat.2011.02.020

  This implementation follows the constant-temperature vacancy-annihilation
  formalism from the above reference.

  For an isothermal hold at temperature T, the vacancy fraction y(t) admits an
  analytic closed form:

    y(t) = y_eq(T) * ( y0 / y_eq(T) )^( exp( -K(T) * t ) )

  where
    - y0        : vacancy fraction at the start of the isothermal segment
                  (e.g., quenched-in value from a prio6r high-T state, or any
                  user-specified initial condition),
    - y_eq(T)   : equilibrium vacancy fraction at temperature T (the long-time limit).

  The relaxation rate is

    K(T) = ( 2 * pi * rho * Dv(T) ) / ( np * f )

  with
    - Dv(T) : vacancy diffusivity at temperature T,
    - rho   : dislocation density,
    - np*d  : jog spacing (np times interplanar spacing d),
    - f     : geometric factor (use f = 0.7272 for bcc in this model).
  */

public:
  /**
   * @brief Construct a rate corrector for isothermal vacancy-annihilation kinetics.
   *
   * Uses the model form
   * `y(t) = y_eq(T) * (y0 / y_eq(T))^(exp(-K(T) * t))`
   * with
   * `K(T) = (2 * pi * rho * Dv(T)) / (np * f)`.
   *
   * @param simCv Vacancy fraction in the supercell.
   * @param diffusivity Vacancy diffusivity term used in the relaxation rate.
   * @param startCv Initial vacancy fraction at the start of the isothermal segment (`y0`).
   * @param eqCv Equilibrium vacancy fraction at T.
   * @param np Dimensionless number (`np * d` is jog spacing).
   * @param rho Dislocation density.
   * @param fValue Geometric factor `f` (default `0.7272` for bcc).
   */
  RateCorrector(
      double simCv,
      double diffusivity,
      double startCv,
      double eqCv,
      double np,
      double rho,
      double fValue = 0.7272)
      : simCv_(simCv),
        startCv_(startCv),
        eqCv_(eqCv),
        kValue_((2.0 * numbers::pi * rho * diffusivity) / (np * fValue))
  {
  }

  [[nodiscard]] inline double GetTimeCorrectionFactor(const double time)
  {
    const double physicalCv = GetCorrectVacancyConcentration(time);
    return simCv_ / physicalCv;
  }

private:
  inline double GetCorrectVacancyConcentration(const double time)
  {
    const double exponent = exp(-kValue_ * time);
    const double base = startCv_ / eqCv_;
    const double cv = eqCv_ * pow(base, exponent);

    return cv;
  }

  /** @brief Species concentrations used for weighted averaging. */
  const double simCv_;
  /** @brief Initial vacancy fraction y0 at the segment start. */
  const double startCv_;
  /** @brief Equilibruim vacancy fraction at T. */
  const double eqCv_;
  /** @brief  K(T) = ( 2 * pi * rho * Dv(T) ) / ( np * f ). */
  const double kValue_;
};

#endif // LMC_LMC_PRED_INCLUDE_RATECORRECTOR_HPP_
