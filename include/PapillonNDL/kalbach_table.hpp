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
#ifndef PAPILLON_NDL_KALBACH_TABLE_H
#define PAPILLON_NDL_KALBACH_TABLE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/interpolation.hpp>
#include <algorithm>
#include <cmath>

namespace pndl {

/**
 * @brief Contains the product Angle-Energy distribution for a single
 *        incident energy, using the Kalbach-Mann representation.
 */
class KalbachTable {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  KalbachTable(const ACE& ace, std::size_t i);

  /**
   * @param energy Outgoing energy grid.
   * @param pdf Probability Density Function for outgoing energy grid.
   * @param cdf Cumulative Density Function for outgoing energy grid.
   * @param R R values as a function of outgoing energy.
   * @param A A values as a function of outgoing energy.
   * @param interp Interpolation method (Histogram or LinLin).
   */
  KalbachTable(const std::vector<double>& energy,
               const std::vector<double>& pdf, const std::vector<double>& cdf,
               const std::vector<double>& R, const std::vector<double>& A,
               Interpolation interp);
  ~KalbachTable() = default;

  double sample_energy(double xi) const {
    auto cdf_it = std::lower_bound(cdf_.begin(), cdf_.end(), xi);
    std::size_t l = std::distance(cdf_.begin(), cdf_it);
    if (xi == *cdf_it) return energy_[l];
    l--;

    // Must account for case where pdf_[l] = pdf_[l+1], which means  that
    // the slope is zero, and m=0. This results in nan for the linear alg.
    // To avoid this, must use histogram for that segment.
    if (interp_ == Interpolation::Histogram || pdf_[l] == pdf_[l + 1])
      return histogram_interp_energy(xi, l);

    return linear_interp_energy(xi, l);
  }

  /**
   * @brief Returns the lowest possible outgoing energy in MeV.
   */
  double min_energy() const { return energy_.front(); }

  /**
   *  @brief Returns the highest possible outgoing energy in MeV.
   */
  double max_energy() const { return energy_.back(); }

  /**
   * @brief Evaluates R for a given outgoing energy.
   * @param E Outgoing energy in MeV.
   */
  double R(double E) const {
    if (E <= energy_.front())
      return R_.front();
    else if (E >= energy_.back())
      return R_.back();
    else {
      auto E_it = std::lower_bound(energy_.begin(), energy_.end(), E);
      std::size_t l = std::distance(energy_.begin(), E_it) - 1;

      if (interp_ == Interpolation::Histogram) {
        return Histogram::interpolate(E, energy_[l], R_[l], energy_[l + 1],
                                      R_[l + 1]);
      } else {
        return LinLin::interpolate(E, energy_[l], R_[l], energy_[l + 1],
                                   R_[l + 1]);
      }
    }
  }

  /**
   * @brief Evaluates A for a given outgoing energy.
   * @param E Outgoing energy in MeV.
   */
  double A(double E) const {
    if (E <= energy_.front())
      return A_.front();
    else if (E >= energy_.back())
      return A_.back();
    else {
      auto E_it = std::lower_bound(energy_.begin(), energy_.end(), E);
      std::size_t l = std::distance(energy_.begin(), E_it) - 1;

      if (interp_ == Interpolation::Histogram) {
        return Histogram::interpolate(E, energy_[l], A_[l], energy_[l + 1],
                                      A_[l + 1]);
      } else {
        return LinLin::interpolate(E, energy_[l], A_[l], energy_[l + 1],
                                   A_[l + 1]);
      }
    }
  }

  /**
   * @breif Evaluates the PDF of scattering with angle mu, and any exit energy.
   * @param mu Cosine of scattering angle.
   */
  double angle_pdf(double mu) const {
    auto mu_p = [](double a, double r, double mu) {
      return 0.5 * (a / std::sinh(a)) *
             (std::cosh(a * mu) + r * std::sinh(a * mu));
    };

    double pdf_out = 0.;

    for (std::size_t i = 0; i < pdf_.size() - 1; i++) {
      if (interp_ == Interpolation::Histogram) {
        pdf_out +=
            mu_p(A_[i], R_[i], mu) * pdf_[i] * (energy_[i + 1] - energy_[i]);
      } else {
        pdf_out += 0.5 * (energy_[i + 1] - energy_[i]) *
                   (mu_p(A_[i], R_[i], mu) * pdf_[i] +
                    mu_p(A_[i + 1], R_[i + 1], mu) * pdf_[i + 1]);
      }
    }

    return pdf_out;
  }

  /**
   * @breif Evaluates the PDF of scattering with angle mu, and exit energy
   *        E_out.
   * @param mu Cosine of scattering angle.
   * @param E_out Exit energy.
   */
  double pdf(double mu, double E_out) const {
    auto mu_p = [](double a, double r, double mu) {
      return 0.5 * (a / std::sinh(a)) *
             (std::cosh(a * mu) + r * std::sinh(a * mu));
    };

    auto E_it = std::lower_bound(energy_.begin(), energy_.end(), E_out);
    if (E_it == energy_.end()) {
      return 0.;
    } else if (E_it == energy_.begin() && E_out < energy_.front()) {
      return 0.;
    }
    std::size_t l = std::distance(energy_.begin(), E_it);
    if (E_out != *E_it) l--;

    if (interp_ == Interpolation::Histogram) {
      double mu_pdf = mu_p(A_[l], R_[l], mu);
      double E_pdf = pdf_[l];
      return mu_pdf * E_pdf;
    } else {
      double f = (E_out - energy_[l]) / (energy_[l + 1] - energy_[l]);
      double out_pdf = 0;

      out_pdf += f * mu_p(A_[l + 1], R_[l + 1], mu) * pdf_[l + 1];
      out_pdf += (1. - f) * mu_p(A_[l], R_[l], mu) * pdf_[l];

      return out_pdf;
    }
  }

  /**
   * @brief Returns a vector of the outgoing energy points.
   */
  const std::vector<double>& energy() const { return energy_; }

  /**
   * @brief Returns a vector for the PDF points corresponding to the
   *        outgoing energy grid.
   */
  const std::vector<double>& pdf() const { return pdf_; }

  /**
   * @brief Returns a vector for the CDF points corresponding to the
   *        outgoing energy grid.
   */
  const std::vector<double>& cdf() const { return cdf_; }

  /**
   * @brief Returns a vector for the values of R corresponding to
   *        the energy grid points.
   */
  const std::vector<double>& R() const { return R_; }

  /**
   * @brief Returns a vector for the values of A corresponding to
   *        the energy grid points.
   */
  const std::vector<double>& A() const { return A_; }

  /**
   * @brief Returns the method of interpolation used for the energy
   *        PDF and CDF, R, and A.
   */
  Interpolation interpolation() const { return interp_; }

  /**
   * @brief Returns the number of outgoing energy points / AngleTables.
   */
  std::size_t size() const { return energy_.size(); }

 private:
  std::vector<double> energy_;
  std::vector<double> pdf_;
  std::vector<double> cdf_;
  std::vector<double> R_;
  std::vector<double> A_;
  Interpolation interp_;

  double histogram_interp_energy(double xi, std::size_t l) const {
    return energy_[l] + ((xi - cdf_[l]) / pdf_[l]);
  }

  double linear_interp_energy(double xi, std::size_t l) const {
    double m = (pdf_[l + 1] - pdf_[l]) / (energy_[l + 1] - energy_[l]);
    return energy_[l] +
           (1. / m) * (std::sqrt(pdf_[l] * pdf_[l] + 2. * m * (xi - cdf_[l])) -
                       pdf_[l]);
  }
};

}  // namespace pndl

#endif
