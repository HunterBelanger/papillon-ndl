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
#ifndef PAPILLON_NDL_LINEARIZE_H
#define PAPILLON_NDL_LINEARIZE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/tabulated_1d.hpp>
#include <functional>
#include <vector>

namespace pndl {

/**
 * @brief Linearizes the function f, keeping the user provided points in
 *        x and y. The resulting linearized function is returnedd as a
 *        Tabulated1D.
 *
 * @warning The values y[i] MUST agree with the values of f(x[i]). This is not
 *          checked however, due to the fact that the function f is allowed to
 *          be discontinuous. A discontinuity is represented with two adjacent
 *          values in the x array which are equal, each having a difference
 *          associated y value. If the values y[i] do not equal f(x[i])
 *          (barring discontinuities), this function may not finish, or may
 *          yield unexpected results.
 *
 * @param x Reference to the vector of x values which will be kept for linear
 *          interpolation. This vector must be sorted.
 * @param y Reference to the vector of y values which will be kept for linear
 *          interpolation. Must be the same length as x.
 * @param f Function to be linearized.
 * @param tolerance Maximum relative absolute error for linear interpolation.
 *                  The default tolerance is 0.001, or 0.1%.
 */
Tabulated1D linearize(const std::vector<double>& x,
                      const std::vector<double>& y,
                      std::function<double(double)> f,
                      double tolerance = 0.001);

/**
 * @brief Linearizes the function f over the interval [x_min, x_max]. The
 *        resulting function is returned as a Tabulated1D.
 *
 * @warning This function may not finish, or may yield unexpected results
 *          if the function f is discontinuous over the interval.
 *
 * @param x_min Minimum x value for linearization.
 * @param x_max Maximum x value for linearization.
 * @param f Function to linearize.
 * @param tolerance Maximum relative absolute error for linear interpolation.
 *                  The default tolerance is 0.001, or 0.1%.
 */
Tabulated1D linearize(double x_min, double x_max,
                      std::function<double(double)> f,
                      double tolerance = 0.001);

}  // namespace pndl

#endif
