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
#ifndef PAPILLON_NDL_DIFFERENCE_1D_H
#define PAPILLON_NDL_DIFFERENCE_1D_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/function_1d.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <memory>

namespace pndl {

/**
 * @brief Class to represent a function which is the difference of two
 *        other functions.
 */
class Difference1D : public Function1D {
 public:
  /**
   * @brief Function will be evaluated as term1(x) - term2(x).
   * @param term1 Pointer to the function for the first term.
   * @param term2 Pointer to the function for the second term.
   */
  Difference1D(std::shared_ptr<Function1D> term1,
               std::shared_ptr<Function1D> term2)
      : term_1_(term1), term_2_(term2) {
    if (!term_1_) {
      std::string mssg = "Term 1 is nullptr.";
      throw PNDLException(mssg);
    }

    if (!term_2_) {
      std::string mssg = "Term 2 is nullptr.";
      throw PNDLException(mssg);
    }
  }

  double operator()(double x) const override final {
    return (*term_1_)(x) - (*term_2_)(x);
  }

  double integrate(double x_low, double x_hi) const override final {
    return term_1_->integrate(x_low, x_hi) - term_2_->integrate(x_low, x_hi);
  }

  /**
   * @brief Returns the first function in the difference.
   */
  const Function1D& term_1() const { return *term_1_; }

  /**
   * @brief Returns the second function in the difference.
   */
  const Function1D& term_2() const { return *term_2_; }

 private:
  std::shared_ptr<Function1D> term_1_;
  std::shared_ptr<Function1D> term_2_;
};

}  // namespace pndl

#endif
