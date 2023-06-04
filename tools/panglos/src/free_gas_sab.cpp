/*
 * Papillon Nuclear Data Library
 * Copyright 2021-2023, Hunter Belanger
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

/**
 * @file
 * @author Hunter Belanger
 */

#include "free_gas_sab.hpp"

#include <cmath>

#include "gauss_kronrod.hpp"

double FreeGasSab::operator()(double a, double b) const {
  const double numerator = std::exp(-((a - b) * (a - b) / (4. * a)) - 0.5 * b);
  const double denominator = std::sqrt(4. * PI * a);
  return numerator / denominator;
}

double FreeGasSab::indefinite_integral_alpha(double a, double b) const {
  b = std::abs(b);

  double erf1 = -1.;
  double erf2 = 1.;

  if (a != 0.) {
    const double inv_2_sqrt_a = 1. / (2. * std::sqrt(a));
    const double erf_arg_1 = (a - b) * inv_2_sqrt_a;
    const double erf_arg_2 = (a + b) * inv_2_sqrt_a;
    erf1 = std::erf(erf_arg_1);
    erf2 = std::erf(erf_arg_2);
  }

  const double exp_b = std::exp(b);
  const double exp_b_2 = std::exp(0.5 * b);

  return 0.5 * exp_b_2 * (-exp_b + erf1 + (exp_b * erf2) + 1.);
}

double FreeGasSab::integrate_alpha(double a_low, double a_hi, double b) const {
  return indefinite_integral_alpha(a_hi, b) -
         indefinite_integral_alpha(a_low, b);
}

double FreeGasSab::integrate_exp_beta(double E, double b_low,
                                      double b_hi) const {
  auto expS = [&, E](double b) {
    return std::exp(-0.5 * b) * this->integrate_alpha(this->min_alpha(E, b),
                                                      this->max_alpha(E, b), b);
  };

  bool flipped = false;
  if (b_low > b_hi) {
    flipped = true;
    double tmp = b_low;
    b_low = b_hi;
    b_hi = tmp;
  }

  GaussKronrodQuadrature<21> GK;
  double integral = GK.integrate(expS, b_low, b_hi, 1.49E-8, 10).first;

  if (flipped) integral = -integral;

  return integral;
}
