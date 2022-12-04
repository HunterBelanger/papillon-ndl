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
#ifndef PAPILLON_NDL_PCTABLE_H
#define PAPILLON_NDL_PCTABLE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/interpolation.hpp>
#include <cmath>
#include <vector>

namespace pndl {

/**
 * @brief Contains a tabulated PDF and CDF for any quantity.
 */
class PCTable {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   * @param normalization Value by which all y-grid points will be divided.
   *                      Default value is 1.
   */
  PCTable(const ACE& ace, std::size_t i, double normalization = 1.);

  /**
   * @param values Vector of values for which the PDF and CDF are provided.
   * @param pdf The Probability Density Function for the provided values.
   * @param cdf The Cumulative Density Function for the provided values.
   * @param interp Interpolation rule for the data. May be either
   *               Histogram or LinLin.
   */
  PCTable(const std::vector<double>& values, const std::vector<double>& pdf,
          const std::vector<double>& cdf, Interpolation interp);
  ~PCTable() = default;

  /**
   * @brief Samples a value from the distribution.
   * @param xi Random value on the interval [0,1).
   */
  double sample_value(double xi) const {
    auto cdf_it = std::lower_bound(cdf_.begin(), cdf_.end(), xi);
    std::size_t l =
        static_cast<std::size_t>(std::distance(cdf_.begin(), cdf_it));
    if (xi == *cdf_it) return values_[l];

    l--;

    // Must account for case where pdf_[l] = pdf_[l+1], which means  that
    // the slope is zero, and m=0. This results in nan for the linear alg.
    // To avoid this, must use histogram for that segment.
    if (interp_ == Interpolation::Histogram || pdf_[l] == pdf_[l + 1])
      return histogram_interp(xi, l);

    return linear_interp(xi, l);
  }

  /**
   * @brief Returns the value of the PDF for the queried value.
   * @param value Value at which to evaluate the PDF.
   */
  double pdf(double value) const {
    if (value < min_value() || value > max_value()) return 0.;

    auto val_it = std::lower_bound(values_.begin(), values_.end(), value);
    std::size_t l =
        static_cast<std::size_t>(std::distance(values_.begin(), val_it));
    if (value == *val_it) return pdf_[l];

    l--;

    if (interp_ == Interpolation::Histogram) return pdf_[l];

    return LinLin::interpolate(value, values_[l], pdf_[l], values_[l + 1],
                               pdf_[l + 1]);
  }

  /**
   * @brief Returns the lowest possible value that can be sampled.
   */
  double min_value() const { return values_.front(); }

  /**
   * @brief Returns the highest possible value that can be sampled.
   */
  double max_value() const { return values_.back(); }

  /**
   * @brief Returns the number of grid points.
   */
  std::size_t size() const { return values_.size(); }

  /**
   * @brief Returns a vector of the value grid points.
   */
  const std::vector<double>& values() const { return values_; }

  /**
   * @brief Returns a vector of the PDF grid points.
   */
  const std::vector<double>& pdf() const { return pdf_; }

  /**
   * @brief Returns a vector of the CDF grid points.
   */
  const std::vector<double>& cdf() const { return cdf_; }

  /**
   * @brief Returns the method of interpolation used for the distribution.
   */
  Interpolation interpolation() const { return interp_; }

 private:
  std::vector<double> values_;
  std::vector<double> pdf_;
  std::vector<double> cdf_;
  Interpolation interp_;

  double histogram_interp(double xi, std::size_t l) const {
    return values_[l] + ((xi - cdf_[l]) / pdf_[l]);
  }

  double linear_interp(double xi, std::size_t l) const {
    double m = (pdf_[l + 1] - pdf_[l]) / (values_[l + 1] - values_[l]);
    return values_[l] +
           (1. / m) * (std::sqrt(pdf_[l] * pdf_[l] + 2. * m * (xi - cdf_[l])) -
                       pdf_[l]);
  }
};

}  // namespace pndl

#endif
