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
#ifndef PAPILLON_NDL_ENERGY_ANGLE_TABLE_H
#define PAPILLON_NDL_ENERGY_ANGLE_TABLE_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/pctable.hpp>
#include <functional>

namespace pndl {

/**
 * @brief Contains the product Angle-Energy distribution for a single
 *        incident energy.
 */
class EnergyAngleTable {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   * @param JED Relative index for finding angular distributions.
   */
  EnergyAngleTable(const ACE& ace, std::size_t i, std::size_t JED);

  /**
   * @param outgoing_energy Outgoing energy grid.
   * @param pdf Probability Density Function for outgoing energy.
   * @param cdf Cumulative Density Function for outgoing energy.
   * @param angle_tables A vector a PCTable, one for each outgoing energy,
   *                     each describing the cosine of the scattering angle.
   * @param interp Interpolation used for the PDF (Histogram or LinLin).
   */
  EnergyAngleTable(const std::vector<double>& outgoing_energy,
                   const std::vector<double>& pdf,
                   const std::vector<double>& cdf,
                   const std::vector<PCTable>& angle_tables,
                   Interpolation interp);

  /**
   * @param outgoing_energy A PCTable cotaining the secondary energy
   *                        distribution.
   * @param angle_tables A vector a PCTable, one for each outgoing energy,
   *                     each describing the cosine of the scattering angle.
   */
  EnergyAngleTable(const PCTable& outgoing_energy,
                   const std::vector<PCTable>& angle_tables);
  ~EnergyAngleTable() = default;

  AngleEnergyPacket sample_angle_energy(
      const std::function<double()>& rng) const {
    double E_out, mu;
    double xi = rng();
    auto cdf_it = std::lower_bound(cdf_.begin(), cdf_.end(), xi);
    std::size_t l = static_cast<std::size_t>(std::distance(cdf_.begin(), cdf_it) - 1);

    // Must account for case where pdf_[l] = pdf_[l+1], which means  that
    // the slope is zero, and m=0. This results in nan for the linear alg.
    // To avoid this, must use histogram for that segment.
    if (interp_ == Interpolation::Histogram || pdf_[l] == pdf_[l + 1]) {
      E_out = histogram_interp_energy(xi, l);
      mu = angles_[l].sample_value(rng());
      if (std::abs(mu) > 1.) mu = std::copysign(1., mu);
      return {mu, E_out};
    }

    E_out = linear_interp_energy(xi, l);

    double f = (xi - cdf_[l]) / (cdf_[l + 1] - cdf_[l]);
    if (f < 0.5)
      mu = angles_[l].sample_value(rng());
    else
      mu = angles_[l + 1].sample_value(rng());

    if (std::abs(mu) > 1.) mu = std::copysign(1., mu);

    return {mu, E_out};
  }

  /**
   * @breif Evaluates the PDF of scattering with angle mu, and any exit energy.
   * @param mu Cosine of scattering angle.
   */
  double angle_pdf(double mu) const {
    double pdf_out = 0.;

    for (std::size_t i = 0; i < pdf_.size() - 1; i++) {
      if (interp_ == Interpolation::Histogram) {
        pdf_out += angles_[i].pdf(mu) * pdf_[i] * (energy_[i + 1] - energy_[i]);
      } else {
        pdf_out += 0.5 * (energy_[i + 1] - energy_[i]) *
                   (angles_[i].pdf(mu) * pdf_[i] +
                    angles_[i + 1].pdf(mu) * pdf_[i + 1]);
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
    auto E_it = std::lower_bound(energy_.begin(), energy_.end(), E_out);
    if (E_it == energy_.end()) {
      return 0.;
    } else if (E_it == energy_.begin() && E_out < energy_.front()) {
      return 0.;
    }
    std::size_t l = static_cast<std::size_t>(std::distance(energy_.begin(), E_it));
    if (E_out != *E_it) l--;

    if (interp_ == Interpolation::Histogram) {
      double mu_pdf = angles_[l].pdf(mu);
      double E_pdf = pdf_[l];
      return mu_pdf * E_pdf;
    } else {
      double f = (E_out - energy_[l]) / (energy_[l + 1] - energy_[l]);
      double out_pdf = 0;

      out_pdf += f * angles_[l + 1].pdf(mu) * pdf_[l + 1];
      out_pdf += (1. - f) * angles_[l].pdf(mu) * pdf_[l];

      return out_pdf;
    }
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
   * @brief Returns the method of interpolation used for the energy
   *        PDF and CDF.
   */
  Interpolation interpolation() const { return interp_; }

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
   * @brief Returns the ith AngleTable which contains the angular
   *        distribution for the ith outgoing energy.
   * @param i Index to the outgoing energy grid.
   */
  const PCTable& angle_table(std::size_t i) const { return angles_[i]; }

  /**
   * @brief Returns the number of outgoing energy points / AngleTables.
   */
  std::size_t size() const { return energy_.size(); }

 private:
  std::vector<double> energy_;
  std::vector<double> pdf_;
  std::vector<double> cdf_;
  std::vector<PCTable> angles_;
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
