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
#ifndef PAPILLON_NDL_TABULATED_1D_H
#define PAPILLON_NDL_TABULATED_1D_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/function_1d.hpp>
#include <PapillonNDL/interpolation.hpp>
#include <vector>

namespace pndl {

/**
 * @brief Interface to represent functions of a single variable which
 *        are represented by a tabulation (TAB1 in ENDF).
 */
class Tabulated1D : public Function1D {
 public:
  virtual ~Tabulated1D() = default;

  /**
   * @brief Returns a vector of the locations in the grid where the
   *        interpolation method changes.
   */
  virtual std::vector<uint32_t> breakpoints() const = 0;

  /**
   * @brief Returns a vector of the interpolation methods for each
   *        segment of the grid.
   */
  virtual std::vector<Interpolation> interpolation() const = 0;

  /**
   * @brief Returns a vector of all x points.
   */
  virtual std::vector<double> x() const = 0;

  /**
   * @brief Returns a vector of all y points.
   */
  virtual std::vector<double> y() const = 0;
};

}  // namespace pndl

#endif
