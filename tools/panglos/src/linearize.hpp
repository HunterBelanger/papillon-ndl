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
#ifndef PANGLOS_LINEARIZE_H
#define PANGLOS_LINEARIZE_H

#include <algorithm>
#include <cmath>
#include <exception>
#include <forward_list>
#include <functional>
#include <stdexcept>
#include <utility>
#include <vector>

struct LinearizedFunction {
  std::vector<double> x, y;
};

inline LinearizedFunction linearize(const std::vector<double>& i_x,
                                    const std::vector<double>& i_y,
                                    const std::function<double(double)>& f,
                                    const double max_rel_dif = 0.01,
                                    const double max_abs_dif = 1.E-7,
                                    const double max_x_abs_dif = 1.E-10) {
  // Do checks on vectors
  if (i_x.size() != i_y.size()) {
    throw std::runtime_error("x and y must have the same length.");
  }

  if (std::is_sorted(i_x.begin(), i_x.end()) == false) {
    throw std::runtime_error("x must be sorted.");
  }

  std::forward_list<double> x(i_x.begin(), i_x.end());
  std::forward_list<double> y(i_y.begin(), i_y.end());

  // Bisect intervals until we are linearly interpolable
  auto xi = x.begin();
  auto xi1 = xi;
  xi1++;
  auto yi = y.begin();
  auto yi1 = yi;
  yi1++;
  while (xi1 != x.end()) {
    // Get mid-point x value. If x[i] == x[i+1], this is a discontinuity, so we
    // continue.
    if (*xi == *xi1 || (*xi1 - *xi) < max_x_abs_dif) {
      xi = xi1++;
      yi = yi1++;
      continue;
    }
    double x_mid = 0.5 * (*xi + *xi1);

    // Get interpolated and real function value
    const double f_interp = 0.5 * (*yi + *yi1);
    const double f_real = f(x_mid);

    // Check tolerance
    const double abs_diff = std::abs(f_interp - f_real);
    const double rel_diff = std::abs(abs_diff / f_real);
    if (rel_diff > max_rel_dif && abs_diff > max_abs_dif && (*xi != x_mid) &&
        (*xi1 != x_mid)) {
      // Add the mid-point
      x.insert_after(xi, x_mid);
      y.insert_after(yi, f_real);
      xi1 = xi;
      xi1++;
      yi1 = yi;
      yi1++;
    } else {
      xi = xi1++;
      yi = yi1++;
    }
  }

  return {{x.begin(), x.end()}, {y.begin(), y.end()}};
}

inline LinearizedFunction linearize(double x_min, double x_max,
                                    const std::function<double(double)>& f,
                                    const double max_rel_dif = 0.01,
                                    const double max_abs_dif = 1.E-7,
                                    const double max_x_abs_dif = 1.E-10) {
  if (x_max <= x_min) {
    throw std::runtime_error("x_max must be larger than x_min.");
  }

  // Initialize vectors for x and y values
  std::vector<double> x, y;
  x.push_back(x_min);
  y.push_back(f(x_min));
  x.push_back(x_max);
  y.push_back(f(x_max));

  return linearize(x, y, f, max_rel_dif, max_abs_dif, max_x_abs_dif);
}

#endif
