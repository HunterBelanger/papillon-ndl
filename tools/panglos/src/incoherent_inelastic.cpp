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

#include "incoherent_inelastic.hpp"

#include <Log.hpp>
#include <boost/hana.hpp>  // Needed for the _c literal for constructing mt4

#include "constants.hpp"
#include "endf_sab.hpp"
#include "interpolator.hpp"
#include "linearize.hpp"
#include "trapezoid.hpp"
using namespace njoy;

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <variant>

static const std::vector<double> njoy_egrid{
    1.E-5,     1.78E-5,   2.5E-5,    3.5E-5,    5.0E-5,   7.0E-5,    1.E-4,
    1.26E-4,   1.6E-4,    2.0E-4,    0.000253,  0.000297, 0.000350,  0.00042,
    0.000506,  0.000615,  0.00075,   0.00087,   0.001012, 0.00123,   0.0015,
    0.0018,    0.00203,   0.002277,  0.0026,    0.003,    0.0035,    0.004048,
    0.0045,    0.005,     0.0056,    0.006325,  0.0072,   0.0081,    0.009108,
    0.01,      0.01063,   0.0115,    0.012397,  0.0133,   0.01417,   0.015,
    0.016192,  0.0182,    0.0199,    0.020493,  0.0215,   0.0228,    0.0253,
    0.028,     0.030613,  0.0338,    0.0365,    0.0395,   0.042757,  0.0465,
    0.050,     0.056925,  0.0625,    0.069,     0.075,    0.081972,  0.09,
    0.096,     0.1035,    0.111573,  0.120,     0.128,    0.1355,    0.145728,
    0.160,     0.172,     0.184437,  0.20,      0.2277,   0.2510392, 0.2705304,
    0.2907501, 0.3011332, 0.3206421, 0.3576813, 0.39,     0.4170351, 0.45,
    0.5032575, 0.56,      0.625,     0.70,      0.78,     0.86,      0.95,
    1.05,      1.16,      1.28,      1.42,      1.55,     1.70,      1.855,
    2.02,      2.18,      2.36,      2.59,      2.855,    3.12,      3.42,
    3.75,      4.07,      4.46,      4.90,      5.35,     5.85,      6.40,
    7.00,      7.65,      8.40,      9.15,      9.85,     10.00};

IncoherentInelastic::IncoherentInelastic(const section::Type<7, 4>& mt4)
    : sab_(), sab_temps_(), awr_(0.), bound_xs_(0.), Emin_(0.), Emax_(0.) {
  auto constants = mt4.constants();

  const int LAT = mt4.LAT();
  const int LASYM = mt4.LASYM();
  const int LLN = constants.LLN();
  awr_ = constants.AWR()[0];

  auto scatteringLaw = mt4.scatteringLaw();

  if (std::holds_alternative<section::Type<7, 4>::TabulatedFunctions>(
          scatteringLaw) == false) {
    throw std::runtime_error(
        "IncoherentInelastic::IncoherentInelastic: No tabulated scattering law "
        "in ENDF file.");
  }

  section::Type<7, 4>::TabulatedFunctions& tsl =
      std::get<section::Type<7, 4>::TabulatedFunctions>(scatteringLaw);

  // Get list of all provided temperatures from the scattering law
  auto lawsForAllTemps = tsl.S();
  sab_temps_ = {lawsForAllTemps[0].T().begin(), lawsForAllTemps[0].T().end()};

  // Get Tab1 for the effective temperature
  section::Type<7, 4>::EffectiveTemperature rawEffectiveTemp =
      mt4.principalEffectiveTemperature();
  Tab1 effectiveTemp =
      makeTab1(rawEffectiveTemp.boundaries(), rawEffectiveTemp.interpolants(),
               rawEffectiveTemp.TMOD(), rawEffectiveTemp.TEFF());

  // Load all S(a,b) scattering laws
  for (std::size_t i = 0; i < sab_temps_.size(); i++) {
    const double T = sab_temps_[i];
    const double Teff = effectiveTemp(T);
    std::unique_ptr<ENDFSab> ENDFST =
        std::make_unique<ENDFSab>(tsl, i, T, Teff, awr_, LAT, LASYM, LLN);
    // std::unique_ptr<TabulatedSab> ST =
    // std::make_unique<LinearizedSab>(*ENDFST); sab_.push_back(std::move(ST));
    sab_.push_back(std::move(ENDFST));
  }

  // Set min and max energy
  Emin_ = 1.E-5;
  Emax_ = std::fmax(5., constants.EMAX());

  // Calculate bound xs for a single nuclide of the principal scatterer
  bound_xs_ = constants.totalFreeCrossSections()[0] * ((awr_ + 1.) / awr_) *
              ((awr_ + 1.) / awr_) / constants.numberAtoms()[0];

  // Write Information
  Log::info("Incoherent Inelastic");
  Log::info("--------------------");
  Log::info("LAT = {}", LAT);
  Log::info("LLN = {}", LLN);
  Log::info("LASYM = {}", LASYM);
  Log::info("AWR = {}", awr_);
  Log::info("Bound XS = {}", bound_xs_);
  Log::info("Num. of Temperatures = {}", sab_temps_.size());
  Log::info("Emax = {} eV", Emax_);
  if (constants.EMAX() < 5.) {
    Log::warning("Evaluation Emax of {} eV was increased to {} eV",
                 constants.EMAX(), Emax_);
  }
  Log::info("");
}

double IncoherentInelastic::ddxs(std::size_t Ti, double Ein, double Eout,
                                 double mu) const {
  if (mu < -1. || mu > 1.) {
    throw std::runtime_error(
        "IncoherentInelastic::ddxs: mu must be in interval [-1, 1].");
  }

  const double T = sab_temps_[Ti];
  const Sab& S = *sab_[Ti];
  const double b = (Eout - Ein) / (KB * T);
  const double a =
      (Eout + Ein - 2. * mu * std::sqrt(Ein * Eout)) / (awr_ * KB * T);
  return (awr_ * bound_xs_ * KB * T / (4. * Ein)) * std::exp(-0.5 * b) *
         S(a, b);
}

double IncoherentInelastic::xs(std::size_t Ti, double Ein) const {
  const double T = sab_temps_[Ti];
  const Sab& S = *sab_[Ti];
  const double b_min = Sab::min_beta(Ein, T);
  const double b_max = Sab::max_beta(Ein, T);
  return (awr_ * bound_xs_ * KB * T / (4. * Ein)) *
         S.integrate_exp_beta(Ein, b_min, b_max);
}

void normalize_pdf_compute_cdf(const std::vector<double>& x,
                               std::vector<double>& pdf,
                               std::vector<double>& cdf) {
  // Make sure x and pdf have same length
  if (x.size() != pdf.size()) {
    throw std::runtime_error("x and pdf must have same size.");
  }

  if (x.size() < 2) {
    throw std::runtime_error("Distribution must have at least 2 points.");
  }

  // Make sure x is sorted
  if (!std::is_sorted(x.begin(), x.end())) {
    throw std::runtime_error("x grid must be sorted.");
  }

  // Make sure PDF is >= 0
  for (std::size_t i = 0; i < pdf.size(); i++) {
    if (pdf[i] < 0.) {
      throw std::runtime_error("pdf must be >= 0.");
    }
  }

  // Initialize the cdf array to be zero
  cdf.resize(pdf.size(), 0.);
  for (auto& c : cdf) c = 0.;

  // Compute the un-normalized CDF
  for (std::size_t i = 1; i < cdf.size(); i++) {
    cdf[i] = cdf[i - 1] + 0.5 * (pdf[i - 1] + pdf[i]) * (x[i] - x[i - 1]);
  }

  // Normalize PDF and CDF
  for (std::size_t i = 0; i < cdf.size(); i++) {
    pdf[i] /= cdf.back();
    cdf[i] /= cdf.back();
  }

  // Make sure CDF starts at 0, and ends at 1, exactly
  cdf.front() = 0.;
  cdf.back() = 1.;
}

void linearize_alpha(const TabulatedSab& S, const double& Ein, const double& b,
                     AlphaDistribution& adist, bool pedantic) {
  // First calculate min and max alpha
  const double alpha_min = S.min_alpha(Ein, b);
  const double alpha_max = S.max_alpha(Ein, b);

  if (alpha_min == alpha_max) {
    // This occurs for beta_min, when Eout = 0. In this case, the idea of a
    // scattering angle is not defined. As such, we just write out the single
    // alpha value, and leave it as so. In practice, one should likely just
    // sample an angle isotropically.
    adist.alpha = {alpha_min};
    adist.pdf = {1.};
    adist.cdf = {1.};
    return;
  }

  std::vector<double> alpha_points;
  alpha_points.push_back(alpha_min);
  for (std::size_t i = 0; i < S.alpha().size(); i++) {
    if (S.alpha()[i] > alpha_min && S.alpha()[i] < alpha_max) {
      alpha_points.push_back(S.alpha()[i]);
    }
  }
  alpha_points.push_back(alpha_max);

  auto apdf = [&S, &b](double a) { return S(a, b); };

  std::vector<double> pdf_points(alpha_points.size(), 0.);
  for (std::size_t i = 0; i < alpha_points.size(); i++) {
    pdf_points[i] = apdf(alpha_points[i]);
  }

  LinearizedFunction alpha_pdf =
      linearize(alpha_points, pdf_points, apdf, 0.01, 1.E-7, 1.E-8);
  adist.alpha = alpha_pdf.x;
  adist.pdf = alpha_pdf.y;

  const double pdf_integral_trap = trapezoid(adist.alpha, adist.pdf);
  if (pedantic) {
    const double pdf_integral = S.integrate_alpha(alpha_min, alpha_max, b);
    const double rel_diff_int =
        std::abs((pdf_integral_trap - pdf_integral) / pdf_integral);
    if (rel_diff_int > 0.05) {
      Log::warning(
          "Alpha distribution for {} eV beta {} has a relative difference of "
          "{} in integral.",
          Ein, b, rel_diff_int);
    }
  }

  // Now we truncate the distribution
  double trunc_integral = 0.;
  const double trunc_int_tol = (1. - ALPHA_INTEGRAL_TOL) * pdf_integral_trap;
  std::size_t i = 1;
  for (i = 1; i < adist.alpha.size(); i++) {
    trunc_integral += 0.5 * (adist.pdf[i] + adist.pdf[i - 1]) *
                      (adist.alpha[i] - adist.alpha[i - 1]);

    if (trunc_integral >= trunc_int_tol && i > 10) {
      break;
    }
  }
  adist.alpha.resize(std::min(i, adist.alpha.size()));
  adist.pdf.resize(adist.alpha.size());
  adist.cdf.resize(adist.alpha.size());

  // Compute the CDF
  normalize_pdf_compute_cdf(adist.alpha, adist.pdf, adist.cdf);

  for (std::size_t i = 0; i < adist.alpha.size(); i++) {
    if (std::isfinite(adist.alpha[i]) == false ||
        std::isfinite(adist.pdf[i]) == false ||
        std::isfinite(adist.cdf[i]) == false) {
      std::stringstream mssg;
      mssg.precision(15);
      mssg << "Bad alpha distribution for Ein = " << Ein << ", b = " << b
           << ".\n";
      mssg << "Alpha = " << adist.alpha[i] << ", pdf = " << adist.pdf[i]
           << ", cdf = " << adist.cdf[i] << ".\n";
      throw std::runtime_error(mssg.str());
    }
  }
}

std::tuple<std::vector<double>, std::vector<double>, double>
determine_initial_beta_grid(const TabulatedSab& S, const double& Ein) {
  // This function gets the points in the beta grid for the outgoing energy
  // distribution. We truncate the beta distribution based on how close we are
  // to the reference integral value, and then normalize this truncated
  // distribution.
  const double beta_min = S.min_beta(Ein);
  const double beta_max = S.max_beta(Ein);
  const double ref_integral = S.integrate_exp_beta(Ein, beta_min, beta_max);

  const auto expS = [&S, &Ein](double b) {
    return std::exp(-0.5 * b) *
           S.integrate_alpha(S.min_alpha(Ein, b), S.max_alpha(Ein, b), b);
  };

  std::vector<double> b, p, c;

  // First add beta_min point
  b.push_back(beta_min);
  p.push_back(expS(beta_min));
  c.push_back(0.);

  if (S.symmetric()) {
    // First add all negative beta points
    for (int ib = static_cast<int>(S.beta().size()) - 1; ib >= 0; ib--) {
      const double beta = -S.beta()[ib];
      if (beta > beta_min && beta < beta_max) {
        c.push_back(c.back() + S.integrate_exp_beta(Ein, b.back(), beta));
        b.push_back(beta);
        p.push_back(expS(beta));
      }
    }
    if (b.back() == 0.) {
      b.pop_back();
      p.pop_back();
      c.pop_back();
    }

    // Now add all positive beta points
    for (std::size_t ib = 0; ib < S.beta().size(); ib++) {
      const double beta = S.beta()[ib];
      if (beta > beta_min && beta < beta_max) {
        c.push_back(c.back() + S.integrate_exp_beta(Ein, b.back(), beta));
        b.push_back(beta);
        p.push_back(expS(beta));

        // Check if we are converged
        if (std::abs(c.back() - ref_integral) <
            BETA_INTEGRAL_TOL * ref_integral) {
          return {b, p, c.back()};
        }
      }
    }
  } else {
    for (std::size_t ib = 0; ib < S.beta().size(); ib++) {
      const double beta = S.beta()[ib];
      if (beta > beta_min && beta < beta_max) {
        c.push_back(c.back() + S.integrate_exp_beta(Ein, b.back(), beta));
        b.push_back(beta);
        p.push_back(expS(beta));

        // Check if we are converged
        if (std::abs(c.back() - ref_integral) <
            BETA_INTEGRAL_TOL * ref_integral) {
          return {b, p, c.back()};
        }
      }
    }
  }

  // We needed all points apparently
  return {b, p, c.back()};
}

void linearize_beta(const TabulatedSab& S, const double& Ein,
                    BetaDistribution& bdist, bool pedantic) {
  // Create grid of beta points which must be in the linearization. This grid
  // comes from the tabulated beta values in the scattering law. We have a
  // special function to do this.
  std::vector<double> beta_points;
  std::vector<double> pdf_points;
  double pdf_integral;
  std::tie(beta_points, pdf_points, pdf_integral) =
      determine_initial_beta_grid(S, Ein);

  auto expS = [&S, &Ein](double b) {
    return std::exp(-0.5 * b) *
           S.integrate_alpha(S.min_alpha(Ein, b), S.max_alpha(Ein, b), b);
  };

  // Perform linearization
  LinearizedFunction beta_pdf =
      linearize(beta_points, pdf_points, expS, 0.01, 1.E-7, 1.E-8);
  bdist.beta = beta_pdf.x;
  bdist.pdf = beta_pdf.y;

  if (pedantic) {
    // After linearization, check that the integration of the PDF using the
    // trapezoid rule is within certain tolerance of the pre-computed integral
    // from when we determined the initial beta grid points.
    const double pdf_integral_trap = trapezoid(bdist.beta, bdist.pdf);
    const double rel_diff_int =
        std::abs((pdf_integral_trap - pdf_integral) / pdf_integral);
    if (rel_diff_int > 0.05) {
      Log::warning(
          "Beta distribution for {} eV has a relative difference of {} in "
          "integral.",
          Ein, rel_diff_int);
    }
  }

  // Compute the CDF
  normalize_pdf_compute_cdf(bdist.beta, bdist.pdf, bdist.cdf);

  for (std::size_t i = 0; i < bdist.beta.size(); i++) {
    if (std::isfinite(bdist.beta[i]) == false ||
        std::isfinite(bdist.pdf[i]) == false ||
        std::isfinite(bdist.cdf[i]) == false) {
      std::stringstream mssg;
      mssg.precision(15);
      mssg << "Bad beta distribution for Ein = " << Ein << ".\n";
      mssg << "Beta = " << bdist.beta[i] << ", pdf = " << bdist.pdf[i]
           << ", cdf = " << bdist.cdf[i] << ".\n";
      throw std::runtime_error(mssg.str());
    }
  }

  // Now, compute the alpha distribution for each beta
  for (std::size_t i = 0; i < bdist.beta.size(); i++) {
    bdist.alpha_dists.emplace_back();
    linearize_alpha(S, Ein, bdist.beta[i], bdist.alpha_dists.back(), pedantic);
  }
}

LinearizedIncoherentInelastic linearize_ii(const IncoherentInelastic& ii,
                                           std::size_t Ti, bool pedantic) {
  LinearizedIncoherentInelastic out;

  // Get the S(a,b) we are evaluating
  const TabulatedSab& S = ii.sab(Ti);

  // First, linearize the II cross section
  Log::info("Linearizing incoherent inelastic cross section...");

  const auto xs =
      std::bind(&IncoherentInelastic::xs, &ii, Ti, std::placeholders::_1);

  std::vector<double> egrid_points;
  egrid_points.push_back(ii.Emin());
  for (std::size_t i = 0; i < njoy_egrid.size(); i++) {
    if (njoy_egrid[i] > ii.Emin() && njoy_egrid[i] < ii.Emax()) {
      egrid_points.push_back(njoy_egrid[i]);
    }
  }
  egrid_points.push_back(ii.Emax());

  std::vector<double> xs_points(egrid_points.size(), 0.);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::size_t i = 0; i < egrid_points.size(); i++) {
    xs_points[i] = xs(egrid_points[i]);
  }

  LinearizedFunction iixs = linearize(egrid_points, xs_points, xs, 0.001);
  out.egrid = iixs.x;
  out.xs = iixs.y;
  Log::info("Number of Energy Grid Points = {}", out.egrid.size());
  Log::info("");

  Log::info("Linearizing incoherent inelastic distribution...");
  // For each incident energy, linearize the beta PDF and all alpha PDF
  out.beta_dists.resize(out.egrid.size());
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
  for (std::size_t ie = 0; ie < out.egrid.size(); ie++) {
    const double Ein = out.egrid[ie];
    linearize_beta(S, Ein, out.beta_dists[ie], pedantic);
  }
  Log::info("Linearization complete.");
  Log::info("");

  return out;
}