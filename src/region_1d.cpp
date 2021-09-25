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
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/region_1d.hpp>
#include <algorithm>

namespace pndl {

Region1D::Region1D(const std::vector<double>& i_x,
                   const std::vector<double>& i_y, Interpolation interp)
    : x_(i_x), y_(i_y), interpolation_(interp), interpolator() {
  if (x_.size() != y_.size()) {
    std::string mssg = "x and y have different sizes. x.size() = " +
                       std::to_string(x_.size()) +
                       " and y.size() = " + std::to_string(y_.size()) + ".";
    throw PNDLException(mssg);
  }

  // Ensure x_ is ordered
  if (!std::is_sorted(x_.begin(), x_.end())) {
    std::string mssg = "x is not sorted.";
    throw PNDLException(mssg);
  }

  // Set interpolator
  if (interpolation_ == Interpolation::Histogram) {
    interpolator = Histogram();
  } else if (interpolation_ == Interpolation::LinLin) {
    interpolator = LinLin();
  } else if (interpolation_ == Interpolation::LinLog) {
    interpolator = LinLog();
  } else if (interpolation_ == Interpolation::LogLin) {
    interpolator = LogLin();
  } else if (interpolation_ == Interpolation::LogLog) {
    interpolator = LogLog();
  }

  auto verify_x_grid = [&i_x](auto& interp) {
    return interp.verify_x_grid(i_x.begin(), i_x.end());
  };
  auto verify_y_grid = [&i_y](auto& interp) {
    return interp.verify_y_grid(i_y.begin(), i_y.end());
  };

  std::visit(verify_x_grid, interpolator);
  std::visit(verify_y_grid, interpolator);
}

bool operator<(const Region1D& R, const double& X) { return R.min_x() < X; }

}  // namespace pndl
