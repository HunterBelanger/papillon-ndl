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
#ifndef PAPILLON_NDL_URRPTABLES_H
#define PAPILLON_NDL_URRPTABLES_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/interpolation.hpp>
#include <PapillonNDL/reaction.hpp>
#include <PapillonNDL/xs_packet.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>
#include <vector>

#include "PapillonNDL/pndl_exception.hpp"

namespace pndl {

/**
 * @brief Class to hold the URR probability tables for a single nuclide,
 *        at a single temperature.
 */
class URRPTables {
 public:
  /**
   * @brief A struct to hold a probability table for a single incident
   *        energy.
   */
  struct PTable {
    /**
     * @brief A struct to hold the cross section values for a single
     *        probability band.
     */
    struct XSBand {
      double total;   /**< Total cross section */
      double elastic; /**< Elastic cross section (MT 2) */
      double fission; /**< Fission cross section (MT 18) */
      double capture; /**< Radiative capture cross section (MT 102) */
      double heating; /**< Heating number */
    };

    std::vector<double> cdf; /**< Probability CDF for cross section bands */
    std::vector<XSBand> xs_bands; /**< Cross section bands */
  };

 public:
  /**
   * @param ace ACE file containing the probability tables.
   * @param total Total cross section of the nuclide.
   * @param disappearance Dissappearance cross section of the nuclide.
   * @param elastic Elastic cross section of the nuclide.
   * @param capture Capture cross section of the nuclide.
   * @param fission Fission cross section of the nuclide.
   * @param heating Heating number CrossSection of the nuclide.
   * @param reactions Vector of all STReaction instances for the nuclide.
   */
  URRPTables(const ACE& ace, const CrossSection& total,
             const CrossSection& disappearance, const CrossSection& elastic,
             const CrossSection& capture, const CrossSection& fission,
             const CrossSection& heating,
             const std::vector<STReaction>& reactions);

  /**
   * @brief Returns true if the PTables are present, and false if not.
   */
  bool is_valid() const { return energy_->size() > 2; }

  /**
   * @brief Calculates the cross section for a given incident energy and
   *        probability. If the incident energy is not within the URR, or if
   *        there are no probability tables, std::nullopt is returned.
   * @param E Incident energy (MeV).
   * @param i Index of E in the global energy grid, for evaluating cross
   *          sections.
   * @param xi Random variable in the interval [0,1).
   */
  std::optional<XSPacket> evaluate_xs(double E, std::size_t i,
                                      double xi) const {
    if (this->is_valid() == false) {
      return std::nullopt;
    }

    if (E < energy_->front() || E > energy_->back()) {
      return std::nullopt;
    }

    // Find the energy index for sampling band
    std::size_t iE = 0;
    double f = 0.;
    auto Eit = std::lower_bound(energy_->begin(), energy_->end(), E);
    if (E == *Eit) {
      if (Eit == --energy_->end()) {
        iE = static_cast<std::size_t>(std::distance(energy_->begin(), Eit) - 1);
        f = 1.;
      } else {
        iE = static_cast<std::size_t>(std::distance(energy_->begin(), Eit));
        f = 0.;
      }
    } else {
      iE = static_cast<std::size_t>(std::distance(energy_->begin(), Eit) - 1);

      if (interp_ == Interpolation::LinLin) {
        f = (E - (*energy_)[iE]) / ((*energy_)[iE + 1] - (*energy_)[iE]);
      } else {
        f = std::log(E / (*energy_)[iE]) /
            std::log((*energy_)[iE + 1] / (*energy_)[iE]);
      }
    }

    // Get reference to the upper and lower PTable
    const PTable& ptable_low = (*ptables_)[iE];
    const PTable& ptable_hi = (*ptables_)[iE + 1];

    // Figure out which band we have sampled
    std::size_t b_low = 0;
    for (b_low = 0; b_low < ptable_low.cdf.size(); b_low++) {
      if (ptable_low.cdf[b_low] >= xi) break;
    }
    std::size_t b_hi = 0;
    for (b_hi = 0; b_hi < ptable_hi.cdf.size(); b_hi++) {
      if (ptable_hi.cdf[b_hi] >= xi) break;
    }

    // XSPacket struct which will contain the returned cross sections
    XSPacket xsout{0., 0., 0., 0., 0., 0., 0.};

    // Evaluate the cross sections depending on interpolation
    const auto& xsb_low = ptable_low.xs_bands[b_low];
    const auto& xsb_hi = ptable_hi.xs_bands[b_hi];
    if (interp_ == Interpolation::LinLin) {
      xsout.elastic = xsb_low.elastic + f * (xsb_hi.elastic - xsb_low.elastic);
      xsout.capture = xsb_low.capture + f * (xsb_hi.capture - xsb_low.capture);
      xsout.fission = xsb_low.fission + f * (xsb_hi.fission - xsb_low.fission);
      xsout.heating = xsb_low.heating + f * (xsb_hi.heating - xsb_low.heating);
    } else {
      if (xsb_low.elastic > 0. && xsb_hi.elastic > 0.) {
        xsout.elastic =
            std::exp(std::log(xsb_low.elastic) +
                     f * std::log(xsb_hi.elastic / xsb_low.elastic));
      } else {
        xsout.elastic = 0.;
      }

      if (xsb_low.capture > 0. && xsb_hi.capture > 0.) {
        xsout.capture =
            std::exp(std::log(xsb_low.capture) +
                     f * std::log(xsb_hi.capture / xsb_low.capture));
      } else {
        xsout.capture = 0.;
      }

      if (xsb_low.fission > 0. && xsb_hi.fission > 0.) {
        xsout.fission =
            std::exp(std::log(xsb_low.fission) +
                     f * std::log(xsb_hi.fission / xsb_low.fission));
      } else {
        xsout.fission = 0.;
      }

      if (xsb_low.heating > 0. && xsb_hi.heating > 0.) {
        xsout.heating =
            std::exp(std::log(xsb_low.heating) +
                     f * std::log(xsb_hi.heating / xsb_low.heating));
      } else {
        xsout.heating = 0.;
      }
    }

    // Check if these are factors. If so, we mulitply by smooth cross sections.
    if (factors_) {
      xsout.elastic *= elastic_(E, i);
      xsout.capture *= capture_(E, i);
      xsout.fission *= fission_(E, i);
      xsout.heating *= heating_(E, i);
    }

    // Set any negatives to zero.
    if (xsout.elastic < 0.) xsout.elastic = 0.;
    if (xsout.capture < 0.) xsout.capture = 0.;
    if (xsout.fission < 0.) xsout.fission = 0.;
    if (xsout.heating < 0.) xsout.heating = 0.;

    // Now get the inelastic portion (if there is any)
    if (inelastic_) {
      xsout.inelastic = inelastic_->evaluate(E, i);
    }

    // Now get other absorption portion (if there is any)
    double other_absorption = 0.;
    if (absorption_) {
      other_absorption = absorption_->evaluate(E, i);
    }

    // Calculate full absorption
    xsout.absorption = xsout.capture + xsout.fission + other_absorption;

    // Calculate total cross section
    xsout.total = xsout.elastic + xsout.inelastic + xsout.absorption;

    return xsout;
  }

  /**
   * @brief Calculates the cross section for a given incident energy and
   *        probability. If the incident energy is not within the URR, or if
   *        there are no probability tables, std::nullopt is returned.
   * @param E Incident energy (MeV).
   * @param xi Random variable in the interval [0,1).
   */
  std::optional<XSPacket> evaluate_xs(double E, double xi) const {
    if (this->is_valid() == false) {
      return std::nullopt;
    }

    // Get the energy index
    std::size_t i = this->elastic_.energy_grid().get_lower_index(E);
    // Call the other method
    return this->evaluate_xs(E, i, xi);
  }

  /**
   * @brief Returns the minimum energy of the URR probability tables.
   */
  double min_energy() const {
    if (energy_->size() == 0) return -1.;
    return energy_->front();
  }

  /**
   * @brief Returns the maximum energy of the URR probability tables.
   */
  double max_energy() const {
    if (energy_->size() == 0) return -1.;
    return energy_->back();
  }

  /**
   * @brief Returns true if provided energy is in the URR energy range.
   * @param E Energy to check, in MeV.
   */
  bool energy_in_range(double E) const {
    if (energy_->size() < 2) return false;
    return (this->min_energy() <= E) && (E <= this->max_energy());
  }

  /**
   * @brief Energies for which a PTable is given.
   */
  const std::vector<double>& energy() const { return *energy_; }

  /**
   * @brief All PTables for the nuclide.
   */
  const std::vector<PTable>& ptables() const { return *ptables_; }

  /**
   * @brief Number of cross section bands in each PTable.
   */
  std::size_t n_xs_bands() const {
    if (ptables_->size() == 0) return 0;
    return ptables_->front().xs_bands.size();
  }

  /**
   * @brief Returns true if the values in the probability tables are factors
   *        which must multiply the smooth cross sections. False is returned
   *        if the actual cross sections are stored.
   */
  bool xs_factors() const { return factors_; }

  /**
   * @brief Returns ture if there the nuclide also has inelastic reactions in
   *        the URR region.
   */
  bool inelastic_competition() const { return inelastic_ != nullptr; }

  /**
   * @brief Returns ture if there the nuclide also has other obsorption
   *        reactions in the URR region.
   */
  bool absorption_competition() const { return absorption_ != nullptr; }

 private:
  Interpolation interp_;
  bool factors_;
  CrossSection total_;  // MT 1
  CrossSection disappearance_;
  CrossSection elastic_;  // MT 2
  CrossSection capture_;  // MT 102
  CrossSection fission_;  // MT 18
  CrossSection heating_;
  std::shared_ptr<CrossSection> inelastic_;
  std::shared_ptr<CrossSection> absorption_;
  std::shared_ptr<std::vector<double>> energy_;
  std::shared_ptr<std::vector<PTable>> ptables_;
};

}  // namespace pndl

#endif
