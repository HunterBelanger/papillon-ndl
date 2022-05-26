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
#include <PapillonNDL/constants.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <span>

namespace pndl {

/**
 * @brief This struct provides methods to facilitate the doppler broadening of
 *        linearly interpolable cross sections, using the SIGMA1 algorithm,
 *        developed by Cullen and Weisbin. Broadening may be reliably performed
 *        regardless of the units of the energy grid. The user may also specify
 *        the desire integration limit, and the approximations to be used at
 *        either end of the energy grid.
 */
struct Sigma1 {
  /**
   * Indicates the units for cross section energy grids.
   * */
  enum class EnergyUnits {
    eV, /**< Energy in units of electron-volts. */
    MeV /**< Energy in units of mega electron-volts. */
  };

  /**
   * Indicates the units for temperatures.
   * */
  enum class TemperatureUnits {
    K,  /**< Temperature in units of Kelvins. */
    eV, /**< Temperature in units of electron-volts. */
    MeV /**< Temperature in units of mega elactron-volts. */
  };

  /**
   * Possible approximations for the cross section at energies outside of the
   * energy grid.
   * */
  enum class Extrapolate {
    Zero, /**< The cross section is zero for energies outside of the provided
             energy grid. */
    Constant, /**< The cross section is constant for energies outside of the
                 provided energy grid. */
    OneOverV, /**< The cross section behaves as 1/v for energies outside of the
                 provided energy grid. */
  };

  /**
   * @brief Computes the alpha broadening parameter given the initial
   *        temperature, the nuclide's atomic weight ratio, and the units for
   *        the temperatures, and the cross section's energy grid.
   * @param T1 The temperature of the provided cross section to be broadend.
   * @param T2 The desired temperature for the cross section to be returned.
   * @param awr The Atomic Weight Ratio of the nuclide to which the provided
   *            cross section is associated.
   * @param tu The temperature units for both T1 and T2.
   * @param eu The energy units for the energy grid of the cross section to be
   *           broadend.
   */
  inline static double alpha(double T1, double T2, double awr,
                             TemperatureUnits tu = TemperatureUnits::K,
                             EnergyUnits eu = EnergyUnits::MeV) {
    if (T2 < T1) {
      std::string mssg = "T2 must be greater than T1.";
      throw PNDLException(mssg);
    }

    // Convert T1 and T2 to be in the same untis as eu.
    switch (eu) {
      case EnergyUnits::eV:
        switch (tu) {
          case TemperatureUnits::eV:
            break;

          case TemperatureUnits::K:
            T1 *= KB * 1.E6;
            T2 *= KB * 1.E6;
            break;

          case TemperatureUnits::MeV:
            T1 *= 1.E6;
            T2 *= 1.E6;
            break;
        }
        break;

      case EnergyUnits::MeV:
        switch (tu) {
          case TemperatureUnits::eV:
            T1 *= 1.E-6;
            T2 *= 1.E-6;
            break;

          case TemperatureUnits::K:
            T1 *= KB;
            T2 *= KB;
            break;

          case TemperatureUnits::MeV:
            break;
        }
        break;
    }

    return awr / (T2 - T1);
  }

  /**
   * @brief Takes a linearly interpolable cross section and computes the doppler
   *        broadened cross section at energy E.
   * @param egrid Energy grid for the initial cross section, in units of MeV.
   * @param xs Cross section values corresponding to the energy grid in egrid.
   * @param E Energy in MeV at which to compute the broadened cross section.
   * @param alpha Constant which is dependent on the temperature change and the
   *              atomic weight ratio of the nuclide in question.
   * @param limit The integration limit for the broadening integrals. In speed
   *              space, the integration will be over the bounds of
   *              y - limit <= x <= y + limit, where x and y are the canonical
   *              speed-like quantities, found in Cullen's original paper.
   * @param low_approx Defines the extrapolation to be used for the cross
   *                   section at energies bellow the lowest value in the
   *                   energy grid.
   * @param hi_approx Defines the extrapolation to be used for the cross
   *                  section at energies above the highest value in the
   *                  energy grid.
   */
  inline static double broaden(std::span<const double> egrid,
                               std::span<const double> xs, const double E,
                               const double alpha, const double limit,
                               const Extrapolate low_approx = Extrapolate::OneOverV,
                               const Extrapolate hi_approx = Extrapolate::Constant) {
    if (alpha < 0.) {
      std::string mssg = "alpha must be greater than zero";
      throw PNDLException(mssg);
    }

    if (xs.size() != egrid.size()) {
      std::string mssg = "xs is a different size than egrid.";
      throw PNDLException(mssg);
    }

    // Check to make sure the desired energy is inside the provided energy grid.
    auto E_it = std::lower_bound(egrid.begin(), egrid.end(), E);
    if (E_it == egrid.end() || (E_it == egrid.begin() && *E_it != egrid[0])) {
      // Desired energy is outside of the provided grid
      std::string mssg = "Desired energy is outside of energy grid.";
      throw PNDLException(mssg);
    }

    // Pre-claulcated needed parameters for the doppler broadening.
    double y = std::sqrt(alpha * E);
    double yy = y * y;
    double inv_y = 1. / y;
    double inv_yy = inv_y * inv_y;

    // Initialize cross section to be returned to zero.
    double sig = 0.; 

    // Initialize Fa, Fb, and H arrays used in the broadening.
    std::array<double, 5> Fa, Fb, H;
    Fa.fill(0.);
    Fb.fill(0.);
    H.fill(0.);

    // Lambda functions to fill Fa, Fb, and H.
    auto fill_f = [](std::array<double, 5>& f, double a) {
      f[0] = 0.5 * std::erfc(a);
      f[1] = (0.5 / std::sqrt(PI)) * std::exp(-a * a);
      f[2] = 0.5 * f[0] + a * f[1];
      f[3] = f[1] + a * a * f[1];
      f[4] = (3. / 2.) * f[2] + a * a * a * f[1];
    };

    auto fill_h = [](const std::array<double, 5>& fa,
                     const std::array<double, 5>& fb,
                     std::array<double, 5>& h) {
      h[0] = fa[0] - fb[0];
      h[1] = fa[1] - fb[1];
      h[2] = fa[2] - fb[2];
      h[3] = fa[3] - fb[3];
      h[4] = fa[4] - fb[4];
    };

    // Find the lower and upper integration limits, based on the provided limit.
    std::size_t low = 0;
    if (y - limit > 0.) {
      auto low_it = std::lower_bound(egrid.begin(), egrid.end(), (y - limit)*(y - limit)/alpha);
      if (low_it != egrid.begin()) {
        low = std::distance(egrid.begin(), low_it) - 1;
      }
    }
    std::size_t hi = egrid.size() - 1;
    auto hi_it = std::lower_bound(egrid.begin(), egrid.end(), (y + limit)*(y + limit)/alpha);
    if (hi_it != egrid.end()) {
      hi = std::distance(egrid.begin(), hi_it);
    }

    //============================================================================
    // Positive Integral
    //----------------------------------------------------------------------------
    if (low == 0 && std::sqrt(alpha*egrid[0]) > y - limit &&
        low_approx == Extrapolate::OneOverV) {
      // Extend the xs as 1/v from 0 to x[0]
      const double x_0 = std::sqrt(alpha*egrid[0]);
      fill_f(Fa, -y);
      fill_f(Fb, x_0 - y);
      fill_h(Fa, Fb, H);
      sig += inv_yy * x_0 * xs[0] * (H[1] + y * H[0]);
    } else if (low == 0 && std::sqrt(alpha*egrid[0]) > y - limit &&
               low_approx == Extrapolate::Constant) {
      // Extend the xs as constant from 0 to x[0]
      const double x_0 = std::sqrt(alpha*egrid[0]);
      fill_f(Fa, -y);
      fill_f(Fb, x_0 - y);
      fill_h(Fa, Fb, H);
      sig += inv_yy * xs[0] * (H[2] + 2. * y * H[1] + yy * H[0]);
    }

    // Do all of the points which are in the grid
#ifdef _OPENMP
#pragma omp simd reduction(+ : sig)
#endif
    for (std::size_t k = low; k < hi; k++) {
      const double x_k = std::sqrt(alpha*egrid[k]);
      const double x_k1 = std::sqrt(alpha*egrid[k + 1]);
      const double sk = (xs[k + 1] - xs[k]) / (x_k1 * x_k1 - x_k * x_k);
      fill_f(Fa, x_k - y);
      fill_f(Fb, x_k1 - y);
      fill_h(Fa, Fb, H);
      const double Ak = inv_yy * (H[2] + 2. * y * H[1] + yy * H[0]);
      const double Bk = inv_yy * H[4] + 4. * inv_y * H[3] + 6. * H[2] +
                  4. * y * H[1] + yy * H[0];
      sig += Ak * (xs[k] - sk * x_k * x_k) + sk * Bk;
    }

    if (hi == egrid.size() - 1 && std::sqrt(alpha*egrid[hi]) < y + limit &&
        hi_approx == Extrapolate::Constant) {
      // Extend the xs as constant from x[N] to infinity
      fill_f(Fa, std::sqrt(alpha*egrid[egrid.size() - 1]) - y);
      sig += xs[egrid.size() - 1] * (inv_yy * Fa[2] + 2. * inv_y * Fa[1] + Fa[0]);
    } else if (hi == egrid.size() - 1 && std::sqrt(alpha*egrid[hi]) < y + limit &&
               hi_approx == Extrapolate::OneOverV) {
      // Extend the xs as 1/v from x[N] to infinity
      fill_f(Fa, std::sqrt(alpha*egrid[egrid.size() - 1]) - y);
      sig += std::sqrt(alpha*egrid[egrid.size() - 1]) * xs[egrid.size() - 1] * inv_yy * (Fa[1] + y * Fa[0]);
    }

    //============================================================================
    // Negative Integral
    //----------------------------------------------------------------------------
    y = -y;
    inv_y = -inv_y;
    // Find the upper integration limit, based on the provided limit.
    hi = egrid.size() - 1;
    hi_it = std::lower_bound(egrid.begin(), egrid.end(), limit*limit/alpha);
    if (hi_it != egrid.end()) {
      hi = std::distance(egrid.begin(), hi_it);
    }

    if (low_approx == Extrapolate::OneOverV) {
      // Extend the xs as 1/v from 0 to x[0]
      const double x_0 = std::sqrt(alpha*egrid[0]);
      fill_f(Fa, -y);
      fill_f(Fb, x_0 - y);
      fill_h(Fa, Fb, H);
      sig -= inv_yy * x_0 * xs[0] * (H[1] + y * H[0]);
    } else if (low_approx == Extrapolate::Constant) {
      // Extend the xs as constant from 0 to x[0]
      double x_0 = std::sqrt(alpha*egrid[0]);
      fill_f(Fa, -y);
      fill_f(Fb, x_0 - y);
      fill_h(Fa, Fb, H);
      sig -= inv_yy * xs[0] * (H[2] + 2. * y * H[1] + yy * H[0]);
    }

#ifdef _OPENMP
#pragma omp simd reduction(- : sig)
#endif
    for (std::size_t k = 0; k < hi; k++) {
      const double x_k = std::sqrt(alpha*egrid[k]);
      const double x_k1 = std::sqrt(alpha*egrid[k + 1]);
      fill_f(Fa, x_k - y);
      fill_f(Fb, x_k1 - y);
      fill_h(Fa, Fb, H);
      const double sk = (xs[k + 1] - xs[k]) / (x_k1 * x_k1 - x_k * x_k);
      const double Ak = inv_yy * (H[2] + 2. * y * H[1] + yy * H[0]);
      const double Bk = inv_yy * H[4] + 4. * inv_y * H[3] + 6. * H[2] +
                  4. * y * H[1] + yy * H[0];
      sig -= Ak * (xs[k] - sk * x_k * x_k) + sk * Bk;
    }

    if (hi == egrid.size() - 1 && std::sqrt(alpha*egrid[hi]) < y + limit &&
        hi_approx == Extrapolate::Constant) {
      // Extend the xs as constant from x[N] to infinity
      fill_f(Fa, std::sqrt(alpha*egrid[egrid.size() - 1]) - y);
      sig -= xs[egrid.size() - 1] * (inv_yy * Fa[2] + 2. * inv_y * Fa[1] + Fa[0]);
    } else if (hi == egrid.size() - 1 && std::sqrt(alpha*egrid[hi]) < y + limit &&
               hi_approx == Extrapolate::OneOverV) {
      // Extend the xs as 1/v from x[N] to infinity
      fill_f(Fa, std::sqrt(alpha*egrid[egrid.size() - 1]) - y);
      sig -= std::sqrt(alpha*egrid[egrid.size() - 1]) * xs[egrid.size() - 1] * inv_yy * (Fa[1] + y * Fa[0]);
    }

    return sig;
  }

  /**
   * @brief Takes a linearly interpolable cross section and computes the doppler
   *        broadened cross section at energy E.
   * @param egrid Energy grid for the initial cross section, in units of MeV.
   * @param xs Cross section values corresponding to the energy grid in egrid.
   * @param E Energy in MeV at which to compute the broadened cross section.
   * @param alpha Constant which is dependent on the temperature change and the
   *              atomic weight ratio of the nuclide in question.
   * @param limit The integration limit for the broadening integrals. In speed
   *              space, the integration will be over the bounds of
   *              y - limit <= x <= y + limit, where x and y are the canonical
   *              speed-like quantities, found in Cullen's original paper.
   * @param low_approx Defines the extrapolation to be used for the cross
   *                   section at energies bellow the lowest value in the
   *                   energy grid.
   * @param hi_approx Defines the extrapolation to be used for the cross
   *                  section at energies above the highest value in the
   *                  energy grid.
   */
  inline static void broaden(std::span<const double>& egrid,
                             std::span<const std::span<const double>>& xs,
                             std::span<double>& xs_out, double E, double alpha,
                             double limit,
                             Extrapolate low_approx = Extrapolate::OneOverV,
                             Extrapolate hi_approx = Extrapolate::Constant) {
    if (alpha < 0.) {
      std::string mssg = "alpha must be greater than zero";
      throw PNDLException(mssg);
    }

    if (xs.size() != xs_out.size()) {
      std::string mssg = "xs and xs_out must have the same size.";
      throw PNDLException(mssg);
    }

    // Vector which contains the start index of each cross-section.
    std::vector<std::size_t> start_indx(xs.size(), 0);

    // Check sizes and get starting indices.
    for (std::size_t r = 0; r < xs.size(); r++) {
      if (xs[r].size() > egrid.size()) {
        std::string mssg =
            "xs[" + std::to_string(r) + "] has more points than the egrid.";
        throw PNDLException(mssg);
      }

      start_indx[r] = egrid.size() - xs[r].size();
    }

    // Check to make sure the desired energy is inside the provided energy grid.
    auto E_it = std::lower_bound(egrid.begin(), egrid.end(), E);
    if (E_it == egrid.end() || (E_it == egrid.begin() && *E_it != egrid[0])) {
      // Desired energy is outside of the provided grid
      std::string mssg = "Desired energy is outside of energy grid.";
      throw PNDLException(mssg);
    }

    // Pre-claulcated needed parameters for the doppler broadening.
    double y = std::sqrt(alpha * E);
    double yy = y * y;
    double inv_y = 1. / y;
    double inv_yy = inv_y * inv_y;

    // Initialize Fa, Fb, and H arrays used in the broadening.
    std::array<double, 5> Fa, Fb, H;
    Fa.fill(0.);
    Fb.fill(0.);
    H.fill(0.);
    for (std::size_t ll = 0; ll < xs_out.size(); ll++) xs_out[ll] = 0.;
    std::vector<double> sk(xs.size(), 0.);

    // Lambda functions to fill Fa, Fb, and H.
    auto fill_f = [](std::array<double, 5>& f, double a) {
      f[0] = 0.5 * std::erfc(a);
      f[1] = (0.5 / std::sqrt(PI)) * std::exp(-a * a);
      f[2] = 0.5 * f[0] + a * f[1];
      f[3] = f[1] + a * a * f[1];
      f[4] = (3. / 2.) * f[2] + a * a * a * f[1];
    };

    auto fill_h = [](const std::array<double, 5>& fa,
                     const std::array<double, 5>& fb,
                     std::array<double, 5>& h) {
      h[0] = fa[0] - fb[0];
      h[1] = fa[1] - fb[1];
      h[2] = fa[2] - fb[2];
      h[3] = fa[3] - fb[3];
      h[4] = fa[4] - fb[4];
    };

    // Find the lower and upper integration limits, based on the provided limit.
    std::size_t low = 0;
    if (y - limit > 0.) {
      auto low_it = std::lower_bound(egrid.begin(), egrid.end(), (y - limit)*(y - limit)/alpha);
      if (low_it != egrid.begin()) {
        low = std::distance(egrid.begin(), low_it) - 1;
      }
    }
    std::size_t hi = egrid.size() - 1;
    auto hi_it = std::lower_bound(egrid.begin(), egrid.end(), (y + limit)*(y + limit)/alpha);
    if (hi_it != egrid.end()) {
      hi = std::distance(egrid.begin(), hi_it);
    }

    //============================================================================
    // Positive Integral
    //----------------------------------------------------------------------------
    if (low == 0 && std::sqrt(alpha*egrid[0]) > y - limit &&
        low_approx == Extrapolate::OneOverV) {
      // Extend the xs as 1/v from 0 to x[0]
      const double x_0 = std::sqrt(alpha*egrid[0]);
      fill_f(Fa, -y);
      fill_f(Fb, x_0 - y);
      fill_h(Fa, Fb, H);
      for (std::size_t ll = 0; ll < xs_out.size(); ll++) {
        // Only use lower extrapolation if this cross-section
        // provides points all the way to the lowest energy.
        // If not, this reaction has a threshold, and we use
        // the Zero extrapolation, no matter what.
        if (start_indx[ll] == 0) {
          xs_out[ll] += inv_yy * x_0 * xs[ll][0] * (H[1] + y * H[0]);
        }
      }
    } else if (low == 0 && std::sqrt(alpha*egrid[0]) > y - limit &&
               low_approx == Extrapolate::Constant) {
      // Extend the xs as constant from 0 to x[0]
      const double x_0 = std::sqrt(alpha*egrid[0]);
      fill_f(Fa, -y);
      fill_f(Fb, x_0 - y);
      fill_h(Fa, Fb, H);
      for (std::size_t ll = 0; ll < xs_out.size(); ll++) {
        if (start_indx[ll] == 0) {
          xs_out[ll] += inv_yy * xs[ll][0] * (H[2] + 2. * y * H[1] + yy * H[0]);
        }
      }
    }

    // Do all of the points which are in the grid
    for (std::size_t k = low; k < hi; k++) {
      const double x_k = std::sqrt(alpha*egrid[k]);
      const double x_k1 = std::sqrt(alpha*egrid[k + 1]);
      for (std::size_t ll = 0; ll < sk.size(); ll++) {
        if (start_indx[ll] <= k) {
          sk[ll] =
              (xs[ll][k + 1 - start_indx[ll]] - xs[ll][k - start_indx[ll]]) /
              (x_k1 * x_k1 - x_k * x_k);
        }
      }
      fill_f(Fa, x_k - y);
      fill_f(Fb, x_k1 - y);
      fill_h(Fa, Fb, H);
      const double Ak = inv_yy * (H[2] + 2. * y * H[1] + yy * H[0]);
      const double Bk = inv_yy * H[4] + 4. * inv_y * H[3] + 6. * H[2] +
                  4. * y * H[1] + yy * H[0];
      for (std::size_t ll = 0; ll < xs_out.size(); ll++) {
        if (start_indx[ll] <= k) {
          xs_out[ll] += Ak * (xs[ll][k - start_indx[ll]] - sk[ll] * x_k * x_k) +
                        sk[ll] * Bk;
        }
      }
    }

    if (hi == egrid.size() - 1 && std::sqrt(alpha*egrid[hi]) < y + limit &&
        hi_approx == Extrapolate::Constant) {
      // Extend the xs as constant from x[N] to infinity
      fill_f(Fa, std::sqrt(alpha*egrid[egrid.size() - 1]) - y);
      for (std::size_t ll = 0; ll < xs_out.size(); ll++)
        xs_out[ll] += xs[ll][egrid.size() - 1] *
                      (inv_yy * Fa[2] + 2. * inv_y * Fa[1] + Fa[0]);
    } else if (hi == egrid.size() - 1 && std::sqrt(alpha*egrid[hi]) < y + limit &&
               hi_approx == Extrapolate::OneOverV) {
      // Extend the xs as 1/v from x[N] to infinity
      fill_f(Fa, std::sqrt(alpha*egrid[egrid.size() - 1]) - y);
      for (std::size_t ll = 0; ll < xs_out.size(); ll++)
        xs_out[ll] += std::sqrt(alpha*egrid[egrid.size() - 1]) * xs[ll][egrid.size() - 1] *
                      inv_yy * (Fa[1] + y * Fa[0]);
    }

    //============================================================================
    // Negative Integral
    //----------------------------------------------------------------------------
    y = -y;
    inv_y = -inv_y;
    // Find the upper integration limit, based on the provided limit.
    hi = egrid.size() - 1;
    hi_it = std::lower_bound(egrid.begin(), egrid.end(), limit*limit/alpha);
    if (hi_it != egrid.end()) {
      hi = std::distance(egrid.begin(), hi_it);
    }

    if (low_approx == Extrapolate::OneOverV) {
      // Extend the xs as 1/v from 0 to x[0]
      const double x_0 = std::sqrt(alpha*egrid[0]);
      fill_f(Fa, -y);
      fill_f(Fb, x_0 - y);
      fill_h(Fa, Fb, H);
      for (std::size_t ll = 0; ll < xs_out.size(); ll++) {
        if (start_indx[ll] == 0) {
          xs_out[ll] -= inv_yy * x_0 * xs[ll][0] * (H[1] + y * H[0]);
        }
      }
    } else if (low_approx == Extrapolate::Constant) {
      // Extend the xs as constant from 0 to x[0]
      const double x_0 = std::sqrt(alpha*egrid[0]);
      fill_f(Fa, -y);
      fill_f(Fb, x_0 - y);
      fill_h(Fa, Fb, H);
      for (std::size_t ll = 0; ll < xs_out.size(); ll++) {
        if (start_indx[ll] == 0) {
          xs_out[ll] -= inv_yy * xs[ll][0] * (H[2] + 2. * y * H[1] + yy * H[0]);
        }
      }
    }

    for (std::size_t k = 0; k < hi; k++) {
      const double x_k = std::sqrt(alpha*egrid[k]);
      const double x_k1 = std::sqrt(alpha*egrid[k + 1]);
      for (std::size_t ll = 0; ll < sk.size(); ll++) {
        if (start_indx[ll] <= k) {
          sk[ll] =
              (xs[ll][k + 1 - start_indx[ll]] - xs[ll][k - start_indx[ll]]) /
              (x_k1 * x_k1 - x_k * x_k);
        }
      }
      fill_f(Fa, x_k - y);
      fill_f(Fb, x_k1 - y);
      fill_h(Fa, Fb, H);
      const double Ak = inv_yy * (H[2] + 2. * y * H[1] + yy * H[0]);
      const double Bk = inv_yy * H[4] + 4. * inv_y * H[3] + 6. * H[2] +
                  4. * y * H[1] + yy * H[0];
      for (std::size_t ll = 0; ll < xs_out.size(); ll++) {
        if (start_indx[ll] <= k) {
          xs_out[ll] -= Ak * (xs[ll][k - start_indx[ll]] - sk[ll] * x_k * x_k) +
                        sk[ll] * Bk;
        }
      }
    }

    if (hi == egrid.size() - 1 && std::sqrt(alpha*egrid[hi]) < y + limit &&
        hi_approx == Extrapolate::Constant) {
      // Extend the xs as constant from x[N] to infinity
      fill_f(Fa, std::sqrt(alpha*egrid[egrid.size() - 1]) - y);
      for (std::size_t ll = 0; ll < xs_out.size(); ll++)
        xs_out[ll] -= xs[ll][egrid.size() - 1] *
                      (inv_yy * Fa[2] + 2. * inv_y * Fa[1] + Fa[0]);
    } else if (hi == egrid.size() - 1 && std::sqrt(alpha*egrid[hi]) < y + limit &&
               hi_approx == Extrapolate::OneOverV) {
      // Extend the xs as 1/v from x[N] to infinity
      fill_f(Fa, std::sqrt(alpha*egrid[egrid.size() - 1]) - y);
      for (std::size_t ll = 0; ll < xs_out.size(); ll++)
        xs_out[ll] -= std::sqrt(alpha*egrid[egrid.size() - 1]) * xs[ll][egrid.size() - 1] *
                      inv_yy * (Fa[1] + y * Fa[0]);
    }
  }
};
}  // namespace pndl

/*
 * References
 * ----------
 *
 * [1] D. E. Cullen and C. R. Weisbin, “Exact Doppler Broadening of Tabulated
 *     Cross Sections,” Nucl Sci Eng, vol. 60, no. 3, pp. 1999–229, 2017,
 *     doi: 10.13182/nse76-1.
 * */
