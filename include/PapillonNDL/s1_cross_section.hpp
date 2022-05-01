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
#ifndef PAPILLON_NDL_S1_CROSS_SECTION_H
#define PAPILLON_NDL_S1_CROSS_SECTION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/energy_grid.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <PapillonNDL/sigma1.hpp>
#include <memory>
#include <span>
#include <sstream>

namespace pndl {

/**
 * @brief Contains the linearly interpolable cross section data for
 *        a single MT, and at a single temperature, which can be
 *        doppler broadened to any temperature above the provided
 *        temperature.
 */
class S1CrossSection {
 public:
  /**
   * @param ace ACE file to take the data from.
   * @param i Index in the XSS block where the cross section starts.
   * @param E_grid EnergyGrid associated with the cross section values.
   * @param get_index Flag to indicate wether the cross section values begin
   *                  at i, or if the energy grid index is at i. Default
   *                  value is true.
   */
  S1CrossSection(const ACE& ace, std::size_t i,
                 const EnergyGrid& E_grid, bool get_index = true);

  /**
   * @param xs Vector containing the cross section values.
   * @param E_grid EnergyGrid to use for the cross section.
   * @param index Starting index in the energy grid.
   * @param temperature Temperature in Kelvin for the provided cross section data.
   * @param awr Atomic weight ratio of the nuclide.
   */
  S1CrossSection(const std::vector<double>& xs, const EnergyGrid& E_grid, std::size_t index, double temperature, double awr);

  /**
   * @param xs Value for the cross section at all points in the provided
   *           energy grid.
   * @param E_grid EnergyGrid to use for the cross section.
   */
  S1CrossSection(double xs, const EnergyGrid& E_grid);

  /**
   * @param xs CrossSection to be used for doppler broadening.
   * @param temperature Temperature in Kelvin for the provided CrossSection.
   * @param awr Atomic weight ratio of the nuclide.
   */
  S1CrossSection(const CrossSection& xs, double temperature, double awr);

  ~S1CrossSection() = default;

  /**
   * @brief Returns value of the cross section at index relative to
   *        the associated energy grid.
   * @param i Index from associated energy grid.
   */
  double operator[](std::size_t i) const {
    return xs_[i]; 
  }

  /**
   * @brief Evaluates the cross section at a given energy. Uses
   *        bisection search.
   * @param T Temperature in Kelvin at which to evaluate the cross section.
   * @param E Energy at which to evaluate the cross section.
   */
  double operator()(double T, double E) const {
    if (xs_.size() == 0) {
      if (E < xs_.energy_grid().min_energy()) return 0.;
      return xs_[0];
    }

    const double diff_T = T - min_temp_;
    
    // We can't doppler boraden to lower temperatures. Throw error.
    if (diff_T < 0.) {
      std::stringstream mssg;
      mssg << "Cannot doppler boraden cross section from " << min_temp_;
      mssg << " Kelvin to " << T << " Kelvin.";
      throw PNDLException(mssg.str());
    }
    
    // If energy is outside of grid, set to max energy.
    if (E > xs_.energy_grid().max_energy()) {
      E = xs_.energy_grid().max_energy();
    }

    if (E < xs_.energy(0)) {
      // If we are below threshold, xs is zero.
      return 0.;
    } else if (E >= xs_.energy_grid().urr_min_energy() ||
               diff_T < 1.) {
      // If We are in the URR, don't doppler broaden.
      // Also, if we are within 1 Kelvin, we will skip doppler broadening.
      return xs_(E); 
    } else {
      // We are in regime which must be doppler boradened.
      const double alpha = Sigma1::alpha(min_temp_, T, awr_);
      std::span<const double> egrid(xs_.energy_grid().grid().begin() + xs_.index(),
                                    xs_.size());
      std::span<const double> xs(xs_.xs().begin(), xs_.size());
      Sigma1::Extrapolate low_approx = Sigma1::Extrapolate::OneOverV;
      if (xs_.index() != 0) {
        low_approx = Sigma1::Extrapolate::Zero; 
      }
      return Sigma1::broaden(egrid, xs, E, alpha, 4., low_approx);
    }
  }

  /**
   * @brief Evaluates the cross section at a given energy. Uses
   *        bisection search.
   * @param T Temperature in Kelvin at which to evaluate the cross section.
   * @param E Energy at which to evaluate the cross section.
   */
  double evaluate(double T, double E) const { return this->operator()(T,E); }

  /**
   * @brief Returns index in the energy grid at which the cross section
   *        values begin.
   */
  std::size_t index() const { return xs_.index(); }

  /**
   * @brief Number of points in the cross section.
   */
  std::size_t size() const { return xs_.size(); }

  /**
   * @brief Returns the ith cross section value.
   */
  double xs(std::size_t i) const { return xs_.xs(i); }

  /**
   * @brief Returns the ith energy value, which corresponds with
   *        the ith cross section value.
   */
  double energy(std::size_t i) const { return xs_[i]; }

  /**
   * @brief Returns the cross section values as a vector of floats.
   */
  const std::vector<double>& xs() const { return xs_.xs(); }

  /**
   * @breif Returns a reference to the EnergyGrid object associated with the
   *        cross section.
   */
  const EnergyGrid& energy_grid() const { return xs_.energy_grid(); }

  /**
   * @brief Returns a copy of the energy grid points for the cross section
   *        as a vector of floats.
   */
  std::vector<double> energy() const { return xs_.energy(); }

  /**
   * @brief Returns the minimum temperature in Kelvin.
   */
  double min_temperature() const { return min_temp_; }

  /**
   * @brief Returns the atomic weight ratio.
   */
  double awr() const { return awr_; }

 private:
  CrossSection xs_;
  double min_temp_;
  double awr_;
};

}  // namespace pndl

#endif
