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
#ifndef PAPILLON_NDL_BETA_ALPHA_TABLE_H
#define PAPILLON_NDL_BETA_ALPHA_TABLE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/pctable.hpp>
#include <cmath>
#include <functional>
#include <vector>

namespace pndl {

/**
 * @brief A struct to hold a sampled alpha and beta.
 */
struct AlphaBetaPacket {
  double alpha; /**< Sampled alpha, for momentum transfer. */
  double beta;  /**< Sampled beta, for energy transfer. */
};

/**
 * @brief Contains the product Beta-Alpha distribution for a single
 *        incident energy. This is used with the DirectSab class.
 */
class BetaAlphaTable {
 public:
  /**
   * @param beta Beta grid.
   * @param pdf Probability Density Function for beta.
   * @param cdf Cumulative Density Function for beta.
   * @param angle_tables A vector a PCTable, one for each beta, each describing
   *                     an alpha distribution.
   */
  BetaAlphaTable(const std::vector<double>& beta,
                 const std::vector<double>& pdf, const std::vector<double>& cdf,
                 const std::vector<PCTable>& alpha_tables);

  ~BetaAlphaTable() = default;

  /**
   * @brief Samples a scattering alpha and beta for incoherent inelastic
   *        sacttering. Returns the values as an AlphaBetaPacket.
   * @param rng Random number generator function.
   */
  AlphaBetaPacket sample_alpha_beta(const std::function<double()>& rng) const {
    double xi = rng();
    auto cdf_it = std::lower_bound(cdf_.begin(), cdf_.end(), xi);
    std::size_t l =
        static_cast<std::size_t>(std::distance(cdf_.begin(), cdf_it) - 1);

    // Must account for case where pdf_[l] = pdf_[l+1], which means  that
    // the slope is zero, and m=0. This results in nan for the linear alg.
    // To avoid this, must use histogram for that segment.
    double b = 0.;
    if (pdf_[l] == pdf_[l + 1]) {
      b = histogram_interp_beta(xi, l);
    } else {
      b = linear_interp_beta(xi, l);
    }

    const double f = (xi - cdf_[l]) / (cdf_[l + 1] - cdf_[l]);
    double a = 0.;
    if (f < 0.5)
      a = alphas_[l].sample_value(rng());
    else
      a = alphas_[l + 1].sample_value(rng());

    return {a, b};
  }

  /**
   * @brief Returns the lowest possible beta.
   */
  double min_beta() const { return beta_.front(); }

  /**
   *  @brief Returns the highest possible beta.
   */
  double max_beta() const { return beta_.back(); }

  /**
   * @brief Returns a vector of the beta points.
   */
  const std::vector<double>& beta() const { return beta_; }

  /**
   * @brief Returns a vector for the PDF points corresponding to the
   *        beta grid.
   */
  const std::vector<double>& pdf() const { return pdf_; }

  /**
   * @brief Returns a vector for the CDF points corresponding to the
   *        beta grid.
   */
  const std::vector<double>& cdf() const { return cdf_; }

  /**
   * @brief Returns the ith PCTable which contains the alpha distribution for
   *        the ith beta.
   * @param i Index to the beta grid.
   */
  const PCTable& alpha_table(std::size_t i) const { return alphas_[i]; }

  /**
   * @brief Returns the number of beta points / PCTables.
   */
  std::size_t size() const { return beta_.size(); }

 private:
  std::vector<double> beta_;
  std::vector<double> pdf_;
  std::vector<double> cdf_;
  std::vector<PCTable> alphas_;

  double histogram_interp_beta(double xi, std::size_t l) const {
    return beta_[l] + ((xi - cdf_[l]) / pdf_[l]);
  }

  double linear_interp_beta(double xi, std::size_t l) const {
    const double m = (pdf_[l + 1] - pdf_[l]) / (beta_[l + 1] - beta_[l]);
    const double arg = pdf_[l] * pdf_[l] + 2. * m * (xi - cdf_[l]);

    return beta_[l] + (1. / m) * (std::sqrt(std::max(arg, 0.)) - pdf_[l]);
  }
};

}  // namespace pndl

#endif
