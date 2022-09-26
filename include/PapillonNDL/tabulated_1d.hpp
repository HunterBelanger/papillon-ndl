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
#ifndef PAPILLON_NDL_TABULATED_1D_H
#define PAPILLON_NDL_TABULATED_1D_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/function_1d.hpp>
#include <PapillonNDL/interpolation.hpp>
#include <algorithm>
#include <cstdint>
#include <span>
#include <vector>

namespace pndl {

/**
 * @brief Interface to represent functions of a single variable which
 *        are represented by a tabulation (TAB1 in ENDF).
 */
class Tabulated1D : public Function1D {
 public:
  /**
   * @param NBT Vector of the breakpoint location.
   * @param INT Vector of the interpolations for each segment.
   * @param x   Vector containing the x grid.
   * @param y   Vector contianing the y grid.
   */
  Tabulated1D(const std::vector<uint32_t>& NBT,
              const std::vector<Interpolation>& INT,
              const std::vector<double>& x, const std::vector<double>& y);

  /**
   * @param interp Interpolation method used for all points.
   * @param x Vector of all x points.
   * @param y Vector of all y points.
   */
  Tabulated1D(Interpolation interp, const std::vector<double>& x,
              const std::vector<double>& y);

  ~Tabulated1D() = default;

  double operator()(double x) const override final {
    if (x <= min_x())
      return (regions_.front())(x);
    else if (x >= max_x())
      return (regions_.back())(x);

    // Get region which contains x
    auto region_it = regions_.begin();
    while (region_it->max_x() < x) region_it++;

    return (*region_it)(x);
  }

  double integrate(double x_low, double x_hi) const override final {
    bool inverted = x_low > x_hi;
    if (inverted) {
      double x_low_tmp = x_low;
      x_low = x_hi;
      x_hi = x_low_tmp;
    }

    // Integration may only be carried out over the function's valid domain
    if (x_low <= min_x())
      x_low = min_x();
    else if (x_low >= max_x())
      x_low = max_x();

    if (x_hi >= max_x())
      x_hi = max_x();
    else if (x_hi <= min_x())
      x_hi = min_x();

    // Get region which contains x_low
    auto region = regions_.begin();
    while (region->max_x() < x_low) region++;

    double integral = 0.;
    double x_low_lim = x_low;
    double x_upp_lim = x_hi;
    bool integrating = true;
    while (integrating) {
      if (x_low_lim < region->min_x()) x_low_lim = region->min_x();
      if (x_upp_lim > region->max_x()) x_upp_lim = region->max_x();

      integral += region->integrate(x_low_lim, x_upp_lim);

      if (x_upp_lim == x_hi)
        integrating = false;
      else {
        x_low_lim = x_upp_lim;
        x_upp_lim = x_hi;
        region++;
      }
    }

    if (inverted) integral = -integral;

    return integral;
  }

  /**
   * @brief Returns a vector of the locations in the grid where the
   *        interpolation method changes.
   */
  const std::vector<uint32_t>& breakpoints() const { return breakpoints_; }

  /**
   * @brief Returns a vector of the interpolation methods for each
   *        segment of the grid.
   */
  const std::vector<Interpolation>& interpolation() const {
    return interpolation_;
  }

  /**
   * @brief Returns a vector of all x points.
   */
  const std::vector<double>& x() const { return x_; }

  /**
   * @brief Returns a vector of all y points.
   */
  const std::vector<double>& y() const { return y_; }

  /**
   * @brief Returns the lowest x value.
   */
  double min_x() const { return x_.front(); }

  /**
   * @brief Returns the highest x value.
   */
  double max_x() const { return x_.back(); }

  /**
   * @brief Linearizes the function to be linearly interpolable to within the
   *        given tolerance.
   *
   * @param tolerance Maximum relative absolute error for linear interpolation.
   *                  The default tolerance is 0.001, or 0.1%.
   */
  void linearize(double tolerance = 0.001);

 private:
  class InterpolationRange {
   public:
    InterpolationRange(Interpolation interp, std::span<const double> x,
                       std::span<const double> y);

    // Methods required by Function1D
    double operator()(double x) const {
      if (x <= min_x())
        return y_.front();
      else if (x >= max_x())
        return y_.back();

      // Get bounding x1 < x < x2
      const auto hi_it = std::lower_bound(x_.begin(), x_.end(), x);
      const auto low_it = hi_it - 1;

      std::size_t i = low_it - x_.begin();

      double x1 = *low_it;
      double x2 = *hi_it;
      double y1 = y_[i];
      double y2 = y_[i + 1];

      return interpolator_.interpolate(x, x1, y1, x2, y2);
    }

    double integrate(double x_low, double x_hi) const {
      bool inverted = x_low > x_hi;
      if (inverted) {
        double x_low_tmp = x_low;
        x_low = x_hi;
        x_hi = x_low_tmp;
      }

      // Integration may only be carried out over the function's valid domain
      if (x_low <= min_x())
        x_low = min_x();
      else if (x_low >= max_x())
        x_low = max_x();

      if (x_hi >= max_x())
        x_hi = max_x();
      else if (x_hi <= min_x())
        x_hi = min_x();

      // Get iterator for lower bound of first interval
      auto low_it = std::lower_bound(x_.begin(), x_.end(), x_low);
      if (*low_it > x_low) low_it--;

      double integral = 0.;
      double x_low_lim = x_low;
      double x_upp_lim = x_hi;
      bool integrating = true;
      while (integrating) {
        auto hi_it = low_it;
        hi_it++;

        std::size_t i = low_it - x_.begin();

        double x1 = *low_it;
        double x2 = *hi_it;
        double y1 = y_[i];
        double y2 = y_[i + 1];

        if (x_low_lim < x1) x_low_lim = x1;
        if (x_upp_lim > x2) x_upp_lim = x2;

        integral +=
            interpolator_.integrate(x_low_lim, x_upp_lim, x1, y1, x2, y2);

        if (x_upp_lim == x_hi)
          integrating = false;
        else {
          x_low_lim = x_upp_lim;
          x_upp_lim = x_hi;
          low_it++;
        }
      }

      if (inverted) integral = -integral;

      return integral;
    }

    std::size_t size() const { return x_.size(); }

    double min_x() const { return x_.front(); }

    double max_x() const { return x_.back(); }

   private:
    std::span<const double> x_;
    std::span<const double> y_;
    Interpolator interpolator_;
  };

  std::vector<uint32_t> breakpoints_;
  std::vector<Interpolation> interpolation_;
  std::vector<double> x_;
  std::vector<double> y_;
  std::vector<InterpolationRange> regions_;
};

}  // namespace pndl

#endif
