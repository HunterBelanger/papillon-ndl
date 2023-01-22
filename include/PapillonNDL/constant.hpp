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
#ifndef PAPILLON_NDL_CONSTANT_H
#define PAPILLON_NDL_CONSTANT_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/function_1d.hpp>

namespace pndl {

/**
 * @brief Implements a function which is a constant value everywhere.
 */
class Constant : public Function1D {
 public:
  /**
   * @param value Value of the function.
   */
  Constant(double value) : value_(value) {}
  ~Constant() = default;

  double operator()(double /*x*/) const override final { return value_; }
  double integrate(double x_low, double x_hi) const override final {
    return value_ * (x_hi - x_low);
  }

 private:
  double value_;
};

}  // namespace pndl

#endif
