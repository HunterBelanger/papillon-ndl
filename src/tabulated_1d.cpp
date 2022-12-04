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
#include <PapillonNDL/linearize.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <cstddef>
#include <cstdint>

namespace pndl {

Tabulated1D::Tabulated1D(const std::vector<uint32_t>& NBT,
                         const std::vector<Interpolation>& INT,
                         const std::vector<double>& x,
                         const std::vector<double>& y)
    : breakpoints_(NBT), interpolation_(INT), x_(x), y_(y), regions_() {
  // Ensure NBT and INT are the same length
  if (NBT.size() != INT.size()) {
    std::string mssg = "NBT and INT have different sizes. NBT.size() = " +
                       std::to_string(NBT.size()) +
                       " and INT.size() = " + std::to_string(INT.size()) + ".";
    throw PNDLException(mssg);
  }

  if (x.size() != y.size()) {
    std::string mssg =
        "x and y have different sizes. x.size() = " + std::to_string(x.size()) +
        " and y.size() = " + std::to_string(y.size()) + ".";
    throw PNDLException(mssg);
  }

  if (!std::is_sorted(x.begin(), x.end())) {
    std::string mssg = "x is not sorted.";
    throw PNDLException(mssg);
  }

  // Make 1D regions of all intervals
  std::size_t low = 0;
  std::size_t hi = 0;
  for (std::size_t i = 0; i < breakpoints_.size(); i++) {
    hi = breakpoints_[i];

    try {
      regions_.push_back(InterpolationRange(
                                            interpolation_[i], {x_.begin() + static_cast<std::ptrdiff_t>(low), x_.begin() + static_cast<std::ptrdiff_t>(hi)},
                                            {y_.begin() + static_cast<std::ptrdiff_t>(low), y_.begin() + static_cast<std::ptrdiff_t>(hi)}));
    } catch (PNDLException& error) {
      std::string mssg = "The i = " + std::to_string(i) +
                         " InterpolationRange could not be constructed when "
                         "building Tabulated1D.";
      error.add_to_exception(mssg);
      throw error;
    }

    low = hi - 1;

    // Check for discontinuity at region boundary
    if (low < x.size() - 1) {
      if (x[low] == x[low + 1]) low++;
    }
  }
}

Tabulated1D::Tabulated1D(Interpolation interp, const std::vector<double>& x,
                         const std::vector<double>& y)
    : breakpoints_(), interpolation_(), x_(x), y_(y), regions_() {
  breakpoints_.push_back(static_cast<uint32_t>(x_.size()));
  interpolation_.push_back(interp);

  if (x.size() != y.size()) {
    std::string mssg =
        "x and y have different sizes. x.size() = " + std::to_string(x.size()) +
        " and y.size() = " + std::to_string(y.size()) + ".";
    throw PNDLException(mssg);
  }

  if (!std::is_sorted(x.begin(), x.end())) {
    std::string mssg = "x is not sorted.";
    throw PNDLException(mssg);
  }

  // Make 1 1D region
  const std::size_t hi = breakpoints_[0];
  try {
    regions_.push_back(InterpolationRange(interpolation_[0],
                                          {x_.begin(), x_.begin() + static_cast<std::ptrdiff_t>(hi)},
                                          {y_.begin(), y_.begin() + static_cast<std::ptrdiff_t>(hi)}));
  } catch (PNDLException& error) {
    std::string mssg =
        "The InterpolationRange could not be constructed when building "
        "Tabulated1D.";
    error.add_to_exception(mssg);
    throw error;
  }
}

void Tabulated1D::linearize(double tolerance) {
  // Check if we are already linear
  if (interpolation_.size() == 1 &&
      interpolation_[0] == Interpolation::LinLin) {
    return;
  }

  std::vector<double> x = x_;
  std::vector<double> y = y_;

  try {
    Tabulated1D new_tab = pndl::linearize(x, y, *this, tolerance);
    *this = new_tab;
  } catch (PNDLException& err) {
    std::string mssg = "Could not linearize Tabulated1D.";
    err.add_to_exception(mssg);
    throw err;
  }
}

Tabulated1D::InterpolationRange::InterpolationRange(Interpolation interp,
                                                    std::span<const double> x,
                                                    std::span<const double> y)
    : x_(x), y_(y), interpolator_(interp) {
  if (x_.size() != y_.size()) {
    std::string mssg = "x and y have different sizes. x.size() = " +
                       std::to_string(x_.size()) +
                       " and y.size() = " + std::to_string(y_.size()) + ".";
    throw PNDLException(mssg);
  }

  if (x_.size() == 0) {
    std::string mssg = "x and y both have a size of zero.";
    throw PNDLException(mssg);
  }

  // Ensure x_ is ordered
  if (!std::is_sorted(x_.begin(), x_.end())) {
    std::string mssg = "x is not sorted.";
    throw PNDLException(mssg);
  }

  interpolator_.verify_x_grid(x_.begin(), x_.end());
  interpolator_.verify_y_grid(y_.begin(), y_.end());
}

}  // namespace pndl
