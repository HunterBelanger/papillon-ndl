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
#ifndef PAPILLON_NDL_MULTI_REGION_1D_H
#define PAPILLON_NDL_MULTI_REGION_1D_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/region_1d.hpp>

namespace pndl {

/**
 * @brief Implementation of a Tabulated1D which is comprised of several
 *        interpolation regions.
 */
class MultiRegion1D : public Tabulated1D {
 public:
  /**
   * @param regions Vector of Region1D objects, one for each interpolation
   *                region.
   */
  MultiRegion1D(const std::vector<Region1D>& regions);

  /**
   * @param NBT Vector of the breakpoint location.
   * @param INT Vector of the interpolations for each segment.
   * @param x   Vector containing the x grid.
   * @param y   Vector contianing the y grid.
   */
  MultiRegion1D(const std::vector<uint32_t>& NBT,
                const std::vector<Interpolation>& INT,
                const std::vector<double>& x, const std::vector<double>& y);
  ~MultiRegion1D() = default;

  // Methods required by Function1D
  double operator()(double x) const override final;
  double integrate(double x_low, double x_hi) const override final;

  // Methods required by Tabulated1D
  std::vector<uint32_t> breakpoints() const override final;
  std::vector<Interpolation> interpolation() const override final;
  std::vector<double> x() const override final;
  std::vector<double> y() const override final;

  // Extra methods
  /**
   * @brief Returns reference to one of the interpolation segments,
   *        which are stored as Region1D objects.
   * @param i Index to the interpolation segment.
   */
  const Region1D& operator[](std::size_t i) const { return regions_[i]; }

  /**
   * @brief Returns the number of interpolation regions.
   */
  std::size_t size() const { return regions_.size(); }

  /**
   * @brief Returns the lowest x value.
   */
  double min_x() const { return regions_.front().min_x(); }

  /**
   * @brief Returns the highest x value.
   */
  double max_x() const { return regions_.back().max_x(); }

 private:
  std::vector<Region1D> regions_;
};

}  // namespace pndl

#endif
