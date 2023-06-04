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

/**
 * @file
 * @author Hunter Belanger
 */

#include "endf_sab.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
#include <range/v3/to_container.hpp>  // Allows us to use ranges::to<cont<type>>(range);
#include <sstream>
#include <stdexcept>

#include "constants.hpp"
#include "gauss_kronrod.hpp"
#include "interpolator.hpp"

using ScatteringFunction =
    section::Type<7, 4>::TabulatedFunctions::ScatteringFunction;

ENDFSab::ENDFSab(section::Type<7, 4>::TabulatedFunctions& TSL,
                 std::size_t indx_T, double T, double Teff, double A, int LAT,
                 int LASYM, int LLN)
    : TabulatedSab(T, Teff, A), alpha_bounds_(), alpha_interps_(), data_() {
  // First, set symmetric_ in the TabulatedSab base class
  symmetric_ = (LASYM == 0);

  // Get the grid of beta values
  beta_ = ranges::to<std::vector<double>>(TSL.betas());

  // Get the breakpoints between beta interpolations
  beta_bounds_ = ranges::to<std::vector<long>>(TSL.boundaries());

  // Get the interpolations between beta values
  beta_interps_.reserve(beta_bounds_.size());
  for (const auto& i : TSL.interpolants()) {
    beta_interps_.emplace_back(Interpolation(i));
  }

  // Get all SaT records, one for each beta
  std::vector<ScatteringFunction> scatteringFuncs =
      ranges::to<std::vector<ScatteringFunction>>(TSL.S());

  // Make sure that beta_ and scatteringFuncs have the same size !
  if (beta_.size() != scatteringFuncs.size()) {
    throw std::runtime_error(
        "Beta Grid and ScatteringFunctions of different sizes");
  }

  // The alpha grid is the same for all values of beta and T. We can then
  // get the alpha grid from the first scattering function.
  alpha_ = ranges::to<std::vector<double>>(scatteringFuncs.front().alphas());
  alpha_bounds_ =
      ranges::to<std::vector<long>>(scatteringFuncs.front().boundaries());
  alpha_interps_.reserve(alpha_bounds_.size());
  for (const auto& i : scatteringFuncs.front().interpolants()) {
    alpha_interps_.emplace_back(Interpolation(i));
  }

  // If LAT = 1, then the alpha and beta grids need to be converted to
  // the true temperature, as they have been stored for room temp
  // T = 0.0253 eV.
  if (LAT == 1) {
    const double C = TROOM / T_;
    for (auto& b : beta_) b *= C;
    for (auto& a : alpha_) a *= C;
  }

  // If ln(S) is stored, we need to change the interpolation ruels.
  // ATTENTION ! The way the ENDF manual describes this is wrong !
  // Written bellow are the proper transformations:
  // LinLog(3) -> LogLog(5)
  // LinLin(2) -> LogLin(4)
  // I don't know what to do with other interpolation rules yet.
  // Maybe the others aren't even used for TSLs. For now, I just
  // rais a warning.
  if (LLN == 1) {
    for (std::size_t i = 0; i < alpha_interps_.size(); i++) {
      if (alpha_interps_[i].interpolation() == Interpolation::LinLog) {
        alpha_interps_[i] = Interpolator(Interpolation::LogLog);
      } else if (alpha_interps_[i].interpolation() == Interpolation::LinLin) {
        alpha_interps_[i] = Interpolator(Interpolation::LogLin);
      } else {
        throw std::runtime_error(
            "Can only have interpolation rule of 2 or 3 with LLN = 1.");
      }
    }
  }

  // Allocate the data_ array
  data_.reallocate({beta_.size(), alpha_.size()});
  for (std::size_t b = 0; b < beta_.size(); b++) {
    std::vector<double> S =
        ranges::to<std::vector<double>>(scatteringFuncs[b].S()[indx_T]);

    // Make sure that alpha_ and S have the same size !
    if (alpha_.size() != S.size()) {
      throw std::runtime_error("Alpha Grid and S grid of different sizes");
    }

    // Add values to data_
    for (std::size_t i = 0; i < alpha_.size(); i++) {
      data_(b, i) = S[i];
    }
  }

  if (LLN == 1) {
    for (std::size_t i = 0; i < data_.size(); i++) {
      data_[i] = std::exp(data_[i]);
    }
  }
}

double ENDFSab::operator()(double a, double b) const {
  const double orig_b = b;
  if (symmetric_ && b < 0.) {
    b = -b;
  }

  // Check that a and b are in range. If not, use the short collision time
  // approximation instead.
  if (a > alpha_.back() || b < beta_.front() || b > beta_.back()) {
    return sct_(a, orig_b);
  }

  // Get the index in the beta grid
  auto beta_it = std::lower_bound(beta_.begin(), beta_.end(), b);
  std::size_t beta_indx = 0;
  if (beta_it == beta_.begin()) {
    beta_indx = 1;
  } else if (beta_it == beta_.end()) {
    throw std::runtime_error("Could not find index for beta.");
  } else {
    beta_indx = std::distance(beta_.begin(), beta_it);
  }
  const double beta_low = beta_[beta_indx - 1];
  const double beta_hi = beta_[beta_indx];

  // Get index of beta interpolator
  std::size_t beta_interp_indx = 0;
  for (std::size_t i = 0; i < beta_bounds_.size() - 1; i++) {
    beta_interp_indx = 0;
    if (beta_indx + 1 <= static_cast<std::size_t>(beta_bounds_[i]) &&
        beta_indx + 1 > static_cast<std::size_t>(beta_bounds_[i + 1])) {
      break;
    }
  }
  const auto& beta_interp = beta_interps_[beta_interp_indx];

  // Check if alpha is bellow the lowest tabulated value
  if (a < alpha_.front()) {
    Interpolator alpha_interp(Interpolation::LinLin);
    if (data_(beta_indx, 0) - data_(beta_indx, 1) > 0.) {
      // S is increasing with decreasing alpha. Use LogLog extrapolation.
      alpha_interp = Interpolator(Interpolation::LogLog);
    } else {
      // S is decreasing with decreasing alpha. Use LogLin extrapolation.
      alpha_interp = Interpolator(Interpolation::LogLin);
    }
    const double S_bh = alpha_interp.interpolate(
        a, alpha_[0], data_(beta_indx, 0), alpha_[1], data_(beta_indx, 1));
    const double S_bl =
        alpha_interp.interpolate(a, alpha_[0], data_(beta_indx - 1, 0),
                                 alpha_[1], data_(beta_indx - 1, 1));
    return beta_interp.interpolate(b, beta_low, S_bl, beta_hi, S_bh);
  }

  // Get the index in the alpha grid
  auto alpha_it = std::lower_bound(alpha_.begin(), alpha_.end(), a);
  std::size_t alpha_indx = 0;
  if (alpha_it == alpha_.begin()) {
    alpha_indx = 1;
  } else if (alpha_it == alpha_.end()) {
    throw std::runtime_error("Could not find index for alpha.");
  } else {
    alpha_indx = std::distance(alpha_.begin(), alpha_it);
  }
  const double alpha_low = alpha_[alpha_indx - 1];
  const double alpha_hi = alpha_[alpha_indx];

  // Get the 4 bounding points for the interpolation
  const double S_bl_al = data_(beta_indx - 1, alpha_indx - 1);
  const double S_bl_ah = data_(beta_indx - 1, alpha_indx);
  const double S_bh_al = data_(beta_indx, alpha_indx - 1);
  const double S_bh_ah = data_(beta_indx, alpha_indx);

  // NJOY does this weird thing, there is any one of the 4 gird points is bellow
  // a certain cutoff value, it will use the SCT approximation, even if a
  // non-zero vlaue is provided in the table. This is weird, but is very much
  // necessary to get values which are similar to those of NJOY.
  if (S_bl_al < SCT_CUTOFF || S_bl_ah < SCT_CUTOFF || S_bh_al < SCT_CUTOFF ||
      S_bh_ah < SCT_CUTOFF) {
    return sct_(a, orig_b);
  }

  // Get the alpha interpolator
  std::size_t alpha_interp_indx = 0;
  for (std::size_t i = 0; i < alpha_bounds_.size() - 1; i++) {
    alpha_interp_indx = 0;
    if (alpha_indx + 1 <= static_cast<std::size_t>(alpha_bounds_[i]) &&
        alpha_indx + 1 > static_cast<std::size_t>(alpha_bounds_[i + 1])) {
      break;
    }
  }
  const auto& alpha_interp = alpha_interps_[alpha_interp_indx];

  const double S_bl =
      alpha_interp.interpolate(a, alpha_low, S_bl_al, alpha_hi, S_bl_ah);
  const double S_bh =
      alpha_interp.interpolate(a, alpha_low, S_bh_al, alpha_hi, S_bh_ah);
  const double S = beta_interp.interpolate(b, beta_low, S_bl, beta_hi, S_bh);

  if (std::isfinite(S) == false) {
    std::stringstream mssg;
    mssg.precision(15);
    mssg << "Calculated S = " << S << " for a = " << a << ", b = " << orig_b
         << ".\n";
    throw std::runtime_error(mssg.str());
  }

  return S;
}

double ENDFSab::integrate_alpha(double a_low, double a_hi, double b) const {
  if (a_low == a_hi) return 0.;

  auto S = [&, b](double a) { return (*this)(a, b); };

  bool flipped = false;
  if (a_low > a_hi) {
    flipped = true;
    double tmp = a_low;
    a_low = a_hi;
    a_hi = tmp;
  }

  std::vector<double> alpha_bounds;
  alpha_bounds.reserve(alpha_.size());
  alpha_bounds.emplace_back(a_low);
  for (std::size_t i = 0; i < alpha_.size(); i++) {
    if (a_low < alpha_[i] && alpha_[i] < a_hi) {
      alpha_bounds.emplace_back(alpha_[i]);
    }
  }
  alpha_bounds.emplace_back(a_hi);

  GaussKronrodQuadrature<21> GK;
  // GaussKronrodQuadrature<15> GK;
  double integral = 0.;

  for (std::size_t i = 0; i < alpha_bounds.size() - 1; i++) {
    integral +=
        GK.integrate(S, alpha_bounds[i], alpha_bounds[i + 1], 1.49E-8, 10)
            .first;
    // integral += GK.integrate(S, alpha_bounds[i], alpha_bounds[i+1]).first;
  }

  if (flipped) integral = -integral;

  return integral;
}

double ENDFSab::integrate_exp_beta(double E, double b_low, double b_hi) const {
  if (b_low == b_hi) return 0.;

  auto expS = [&, E](double b) {
    return std::exp(-0.5 * b) * this->integrate_alpha(this->min_alpha(E, b),
                                                      this->max_alpha(E, b), b);
  };

  bool flipped = false;
  if (b_low > b_hi) {
    flipped = true;
    double tmp = b_low;
    b_low = b_hi;
    b_hi = tmp;
  }

  std::vector<double> beta_bounds;
  beta_bounds.reserve(beta_.size());
  beta_bounds.emplace_back(b_low);
  for (std::size_t i = 0; i < beta_.size(); i++) {
    if (b_low < beta_[i] && beta_[i] < b_hi) {
      beta_bounds.emplace_back(beta_[i]);
    }
  }
  beta_bounds.emplace_back(b_hi);

  GaussKronrodQuadrature<21> GK;
  // GaussKronrodQuadrature<15> GK;
  double integral = 0.;

  for (std::size_t i = 0; i < beta_bounds.size() - 1; i++) {
    integral +=
        GK.integrate(expS, beta_bounds[i], beta_bounds[i + 1], 1.49E-8, 10)
            .first;
    // integral += GK.integrate(expS, beta_bounds[i], beta_bounds[i+1]).first;
  }

  if (flipped) integral = -integral;

  return integral;
}
