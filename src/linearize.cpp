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

#include <PapillonNDL/interpolation.hpp>
#include <PapillonNDL/linearize.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>
#include <vector>

namespace pndl {

Tabulated1D linearize(const std::vector<double>& i_x,
                      const std::vector<double>& i_y,
                      std::function<double(double)> f, double tolerance) {
  std::vector<double> x = i_x;
  std::vector<double> y = i_y;

  // Do checks on vectors
  if (x.size() != y.size()) {
    std::string mssg = "x and y must have the same length.";
    throw PNDLException(mssg);
  }

  if (std::is_sorted(x.begin(), x.end()) == false) {
    std::string mssg = "x must be sorted.";
    throw PNDLException(mssg);
  }

  // Bisect intervals until we are linearly interpolable
  std::size_t i = 0;
  while (i < (x.size() - 1)) {
    // Get mid-point x value. If x[i] == x[i+1], this is a discontinuity, so we
    // continue.
    if (std::nextafter(x[i], x[i + 1]) == x[i + 1]) {
      i++;
      continue;
    }
    double x_mid = 0.5 * (x[i] + x[i + 1]);

    // Get interpolated and real function value
    const double f_interp = 0.5 * (y[i] + y[i + 1]);
    const double f_real = f(x_mid);

    // Check tolerance
    const double rel_diff = std::abs((f_interp - f_real) / f_real);
    if (rel_diff > tolerance) {
      // Add the mid-point
      auto xp = x.begin() + static_cast<std::ptrdiff_t>(i) + 1;
      auto yp = y.begin() + static_cast<std::ptrdiff_t>(i) + 1;
      x.insert(xp, x_mid);
      y.insert(yp, f_real);
    } else {
      i++;
    }
  }

  x.shrink_to_fit();
  y.shrink_to_fit();

  return Tabulated1D(Interpolation::LinLin, x, y);
}

Tabulated1D linearize(double x_min, double x_max,
                      std::function<double(double)> f, double tolerance) {
  if (x_max <= x_min) {
    std::string mssg = "x_max must be larger than x_min.";
    throw PNDLException(mssg);
  }

  // Initialize vectors for x and y values
  std::vector<double> x, y;
  x.push_back(x_min);
  y.push_back(f(x_min));
  x.push_back(x_max);
  y.push_back(f(x_max));

  return linearize(x, y, f, tolerance);
}

}  // namespace pndl
