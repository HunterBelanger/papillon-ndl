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
#ifndef PAPILLON_NDL_REGION_1D_H
#define PAPILLON_NDL_REGION_1D_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/interpolation.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <memory>
#include <variant>

namespace pndl {

/**
 * @brief Implementation of a Tabulated1D which has only one interpolation
 *        retion.
 */
class Region1D : public Tabulated1D {
 public:
  /**
   * @param i_x Vector of all x points.
   * @param i_y Vector of all y points.
   * @param interp Interpolation method used for all points.
   */
  Region1D(const std::vector<double>& i_x, const std::vector<double>& i_y,
           Interpolation interp);

  // Methods required by Function1D
  double operator()(double x) const override final {
    if (x <= min_x())
      return y_.front();
    else if (x >= max_x())
      return y_.back();

    // Get bounding x1 < x < x2
    auto low_it = std::lower_bound(x_.begin(), x_.end(), x);
    low_it--;

    auto hi_it = low_it;
    hi_it++;

    std::size_t i = low_it - x_.begin();

    double x1 = *low_it;
    double x2 = *hi_it;
    double y1 = y_[i];
    double y2 = y_[i + 1];

    auto doInterp = [&x, &x1, &x2, &y1, &y2](auto& interp) {
      return interp.interpolate(x, x1, y1, x2, y2);
    };
    return std::visit(doInterp, interpolator);
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

      auto doIntegrl = [&x_low_lim, &x_upp_lim, &x1, &x2, &y1,
                        &y2](auto& interp) {
        return interp.integrate(x_low_lim, x_upp_lim, x1, y1, x2, y2);
      };
      integral += std::visit(doIntegrl, interpolator);

      // integral += interpolator.integrate(x_low_lim, x_upp_lim, x1, y1, x2,
      // y2);

      if (x_upp_lim == x_hi)
        integrating = false;
      else {
        x_low_lim = x_upp_lim;
        x_upp_lim = x_hi;
        low_it++;
      }
    }

    if (inverted) integral *= -1.;

    return integral;
  }

  // Methods required by Tabulated1D
  std::vector<uint32_t> breakpoints() const override final {
    return {static_cast<uint32_t>(x_.size())};
  }
  std::vector<Interpolation> interpolation() const override final {
    return {interpolation_};
  }
  std::vector<double> x() const override final { return x_; }
  std::vector<double> y() const override final { return y_; }

  /**
   * @brief Returns the number of (x,y) pairs.
   */
  std::size_t size() const { return x_.size(); }

  /**
   * @brief Returns the lowest x value.
   */
  double min_x() const { return x_.front(); }

  /**
   * @brief Returns the highest x value.
   */
  double max_x() const { return x_.back(); }

 private:
  std::vector<double> x_;
  std::vector<double> y_;
  Interpolation interpolation_;
  std::variant<Histogram, LinLin, LinLog, LogLin, LogLog> interpolator;
};

// This operator overload is provided only to accomodate the
// std::lower_bound algorithm, in the MultiRegion1D class
bool operator<(const Region1D& R, const double& X);

}  // namespace pndl

#endif
