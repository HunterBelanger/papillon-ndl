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
#ifndef PANGLOS_SAB_H
#define PANGLOS_SAB_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <cmath>
#include <exception>
#include <stdexcept>

#include "constants.hpp"

/**
 * @brief This abstract class acts as an interface for \f$S(\alpha,\beta)\f$
 *        scattering laws. The two dimension-less parameters \f$\alpha\f$ and
 *        \f$\beta\f$ are defined as
 *        \f[ \beta = \frac{E' - E}{k T} \f]
 *        \f[ \alpha = \frac{E' + E - 2\mu\sqrt{EE'}}{A k T} \f]
 *        where \f$E\f$ is the incident energy, \f$E'\f$ is the exit energy,
 *        \f$\mu\f$ is the cosine of the scattering angle, \f$T\f$ is the
 *        temerature of the material, and \f$k\f$ is the Boltzmann constant.
 **/
class Sab {
 public:
  /**
   * @param T Temerature of the scattering law.
   * @param A The atomic weight ratio of the isotope for the scattering law.
   **/
  Sab(double T, double A) : T_(T), A_(A) {
    if (T_ <= 0.) {
      throw std::runtime_error("Temperature for Sab must be > 0.");
    }

    if (A_ <= 0.) {
      throw std::runtime_error("AWR for Sab must be > 0.");
    }
  }

  /**
   * @brief Evaluate the \f$S(\alpha,\beta)\f$ function.
   * @param a Value of \f$\alpha\f$.
   * @param b Value of \f$\beta\f$.
   **/
  virtual double operator()(double a, double b) const = 0;

  /**
   * @brief Evaluates the integral
            \f[
                \int_{\alpha_\text{low}}^{\alpha_\text{hi}}
                     S(\alpha,\beta)
                d\alpha
            \f]
   * @param a_low Lower limit of integration (\f$\alpha_\text{low}\f$).
   * @param a_hi Upper limit of integration (\f$\alpha_\text{hi}\f$).
   * @param b Value of \f$\beta\f$.
   **/
  virtual double integrate_alpha(double a_low, double a_hi, double b) const = 0;

  /**
   * @brief Evaluates the double integral
            \f[
                \int_{\beta_\text{low}(E)}^{\beta_\text{hi}(E)} e^{-\beta/2}
                \int_{\alpha_\text{min}(\beta)}^{\alpha_\text{max}(\beta)}
                     S(\alpha,\beta)
                d\alpha
                d\beta
            \f]
            This integral is typically used when reconstructing the integral
            cross section for incoherent inelastic scattering.
   * @param E Incident energy in eV.
   * @param b_low Lower limit of integration (\f$\beta_\text{low}\f$).
   * @param b_hi Upper limit of integration (\f$\beta_\text{hi}\f$).
   **/
  virtual double integrate_exp_beta(double E, double b_low,
                                    double b_hi) const = 0;

  /**
   * @brief Determines the minimum value of \f$\beta\f$, by setting
   *        \f$E' = 0\f$.
   * @param E Incident energy in eV.
   */
  double min_beta(double E) const { return -E / (KB * T_); }

  /**
   * @brief Determines the minimum value of \f$\beta\f$, by setting
   *        \f$E' = 0\f$.
   * @param E Incident energy in eV.
   * @param T Temerature of the moderator in K.
   */
  static double min_beta(double E, double T) { return -E / (KB * T); }

  /**
   * @brief Determines the maximum value of \f$\beta\f$. There is no theoretical
   *        limit for the maximum energy transfer, and is is hard-coded to 20.
   *        Energy transters larger than this would be exceptionally rare
   *        events.
   * @param E Incident energy in eV.
   */
  double max_beta(double E) const {
    (void)E;
    return 20.;
  }

  /**
   * @brief Determines the maximum value of \f$\beta\f$. There is no theoretical
   *        limit for the maximum energy transfer, and is is hard-coded to 20.
   *        Energy transters larger than this would be exceptionally rare
   *        events.
   * @param E Incident energy in eV.
   * @param T Temerature of the moderator in K.
   */
  static double max_beta(double E, double T) {
    (void)E;
    (void)T;
    return 20.;
  }

  /**
   * @brief Determines the minimum value of \f$\alpha\f$ for a given
   *        \f$\beta\f$, using \f$\mu = 1\f$.
   * @param E Incident energy in eV.
   * @param b The value of \f$\beta\f$.
   */
  double min_alpha(double E, double b) const {
    const double sqrt_num =
        std::sqrt(E) - std::sqrt(std::fmax(E + b * KB * T_, 0.));
    const double num = sqrt_num * sqrt_num;
    const double denom = A_ * KB * T_;
    return num / denom;
  }

  /**
   * @brief Determines the maximum value of \f$\alpha\f$ for a given
   *        \f$\beta\f$, using \f$\mu = -1\f$.
   * @param E Incident energy in eV.
   * @param b The value of \f$\beta\f$.
   */
  double max_alpha(double E, double b) const {
    const double sqrt_num =
        std::sqrt(E) + std::sqrt(std::fmax(E + b * KB * T_, 0.));
    const double num = sqrt_num * sqrt_num;
    const double denom = A_ * KB * T_;
    return num / denom;
  }

  double temperature() const { return T_; }

  double atomic_weight_ratio() const { return A_; }

 protected:
  double T_, A_;
};

#endif
