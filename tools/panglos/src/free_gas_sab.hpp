/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2022, Hunter Belanger
 *
 * hunter.belanger@gmail.com
 *
 * This file is part of the Papillon Nuclear Data Library (PapillonNDL).
 *
 * PapillonNDL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PapillonNDL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PapillonNDL. If not, see <https://www.gnu.org/licenses/>.
 *
 * */
#ifndef PANGLOS_FREE_GAS_SAB_H
#define PANGLOS_FREE_GAS_SAB_H

/**
 * @file
 * @author Hunter Belanger
 */

#include "constants.hpp"
#include "sab.hpp"

/**
   @brief This class is used for the Free Gas Approximation which
          has the following form:
          \f[
               S(\alpha,\beta) =
               \frac{\exp\left(-\frac{(\alpha-\beta)^2}
                                     {4\alpha}
                               -\frac{\beta}{2}
                         \right)
                    }
                    {\sqrt{4\pi\alpha}}
          \f]
**/
class FreeGasSab : public Sab {
 public:
  /**
   * @param T Actual temperature in K.
   * @param A Atomic weigh ratio of the isotope.
   **/
  FreeGasSab(double T, double A)
      : Sab(T, A) {}

  double operator()(double a, double b) const override final;

  double integrate_alpha(double a_low, double a_hi,
                         double b) const override final;

  double integrate_exp_beta(double E, double b_low,
                            double b_hi) const override final;

 private:
  double indefinite_integral_alpha(double a, double b) const;
};

#endif
