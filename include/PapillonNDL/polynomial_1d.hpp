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
#ifndef PAPILLON_NDL_POLYNOMIAL_H
#define PAPILLON_NDL_POLYNOMIAL_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/function_1d.hpp>
#include <vector>

namespace pndl {

/**
 * @brief A function of a single variable, represented by polynomial
 *        coefficients.
 */
class Polynomial1D : public Function1D {
 public:
  /**
   * @param coeffs Vector containing all of the polynomial coefficients.
   */
  Polynomial1D(const std::vector<double>& coeffs);
  ~Polynomial1D() = default;

  double operator()(double x) const override final;
  double integrate(double x_low, double x_hi) const override final;

  /**
   * @brief Returns the order of the polynomial.
   */
  std::size_t order() const { return coefficients_.size() - 1; }

  /**
   * @brief Returns the coefficient for the ith order term.
   * @param i Order of the term.
   */
  double coefficient(std::size_t i) const { return coefficients_[i]; }

 private:
  std::vector<double> coefficients_;
};

}  // namespace pndl

#endif
