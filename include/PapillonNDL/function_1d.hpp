/*
 * Papillon Nuclear Data Library
 * Copyright 2021, Hunter Belanger
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
#ifndef PAPILLON_NDL_FUNCTION_1D_H
#define PAPILLON_NDL_FUNCTION_1D_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <memory>

namespace pndl {

/**
 * @brief Interface to represent functions of a single variable.
 */
class Function1D : public std::enable_shared_from_this<Function1D> {
 public:
  virtual ~Function1D() = default;

  /**
   * @brief Evaluates the function for a given value.
   * @param x Value at which to evaluate the function.
   */
  virtual double operator()(double x) const = 0;

  /**
   * @brief Evaluates the function for a given value.
   * @param x Value at which to evaluate the function.
   */
  double evaluate(double x) const { return (*this)(x); }

  /**
   * @brief Computes the definite integral of the function between
   *        two values.
   * @param x_low Lower bound of integration.
   * @param x_hi Upper bound of integration.
   */
  virtual double integrate(double x_low, double x_hi) const = 0;
};

}  // namespace pndl

#endif
