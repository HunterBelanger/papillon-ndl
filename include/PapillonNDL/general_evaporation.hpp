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
#ifndef PAPILLON_NDL_GENERAL_EVAPORATION_H
#define PAPILLON_NDL_GENERAL_EVAPORATION_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <memory>
#include <optional>

namespace pndl {

/**
 * @brief Energy distribution represented as a general evaporation spectrum.
 */
class GeneralEvaporation : public EnergyLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  GeneralEvaporation(const ACE& ace, std::size_t i);

  /**
   * @param temperature Tabulated function for nuclear temperature.
   * @param bounds Equiprobable bounds for X.
   */
  GeneralEvaporation(std::shared_ptr<Tabulated1D> temperature,
                     const std::vector<double>& bounds);
  ~GeneralEvaporation() = default;

  double sample_energy(double E_in,
                       std::function<double()> rng) const override final;

  std::optional<double> pdf(double E_in, double E_out) const override final;

  /**
   * @brief Returns the table containg the effective temperature
   *        as a function of incoming energy.
   */
  const Tabulated1D& temperature() const { return *temperature_; }

  /**
   * @brief Returns the the bin boundaries in a vector.
   */
  const std::vector<double>& bin_bounds() const { return bin_bounds_; }

 private:
  std::shared_ptr<Tabulated1D> temperature_;
  std::vector<double> bin_bounds_;
};

}  // namespace pndl

#endif
