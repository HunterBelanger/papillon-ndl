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
#include <PapillonNDL/polynomial_1d.hpp>

namespace pndl {

Polynomial1D::Polynomial1D(const std::vector<double>& coeffs)
    : coefficients_(coeffs) {}

double Polynomial1D::operator()(double x) const {
  double value = 0.;
  for (auto coeff = coefficients_.rbegin(); coeff != coefficients_.rend();
       coeff++) {
    value = (value * x) + *coeff;
  }
  return value;
}

double Polynomial1D::integrate(double x_low, double x_hi) const {
  double integral_hi = 0.;
  double integral_low = 0.;

  double exponent = static_cast<double>(order() + 1);
  for (auto coeff_it = coefficients_.rbegin(); coeff_it != coefficients_.rend();
       coeff_it++) {
    double coeff = *coeff_it / exponent;
    integral_hi = (integral_hi * x_hi) + coeff;
    integral_low = (integral_low * x_low) + coeff;
    exponent -= 1.;
  }
  // Do one extra where coeff = 0
  integral_hi = integral_hi * x_hi;
  integral_low = integral_low * x_low;

  return integral_hi - integral_low;
}

}  // namespace pndl
