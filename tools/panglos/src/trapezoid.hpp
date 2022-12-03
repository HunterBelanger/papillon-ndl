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
#ifndef PANGLOS_TRAPEZOID_H
#define PANGLOS_TRAPEZOID_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <exception>
#include <stdexcept>
#include <vector>

inline double trapezoid(const std::vector<double>& x, const std::vector<double>& y) {
  // Do checks on vectors
  if (x.size() != y.size()) {
    throw std::runtime_error("x and y must have the same length.");
  }

  if (x.size() < 2) {
    throw std::runtime_error("x and y must have at least two elements.");
  }

  if (std::is_sorted(x.begin(), x.end()) == false) {
    throw std::runtime_error("x must be sorted.");
  }

  double integral = 0;

  for (std::size_t i = 0; i < x.size()-1; i++) {
    integral += 0.5*(y[i] + y[i+1]) * (x[i+1] - x[i]);
  }

  return integral;
}

#endif
