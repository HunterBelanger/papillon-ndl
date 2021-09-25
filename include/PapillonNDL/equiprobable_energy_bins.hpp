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
#ifndef PAPILLON_NDL_EQUIPROBABLE_ENERGY_BINS_H
#define PAPILLON_NDL_EQUIPROBABLE_ENERGY_BINS_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <PapillonNDL/energy_law.hpp>
#include <optional>

namespace pndl {

/**
 * @brief Energy distribution which is provided as equiprobable energy bins.
 */
class EquiprobableEnergyBins : public EnergyLaw {
 public:
  /**
   * @param ace ACE file to take data from.
   * @param i Starting index of distribution in the XSS array.
   */
  EquiprobableEnergyBins(const ACE& ace, std::size_t i);

  /**
   * @param incoming_energy Incoming energy grid.
   * @param bin_bounds A vector of vectors containing bin bounds, one for
   *                   each point in the incoming energy grid.
   */
  EquiprobableEnergyBins(const std::vector<double>& incoming_energy,
                         const std::vector<std::vector<double>>& bin_bounds);
  ~EquiprobableEnergyBins() = default;

  double sample_energy(double E_in,
                       std::function<double()> rng) const override final;

  std::optional<double> pdf(double E_in, double E_out) const override final;

  /**
   * @brief Returns a vector of the grid of incoming energy points for which
   *        an equiprobable bin set is stored.
   */
  const std::vector<double>& incoming_energy() const {
    return incoming_energy_;
  }

  /**
   * @brief Returns the ith set of bin boundaries as a vector.
   * @param i Index for the incoming energy grid.
   */
  const std::vector<double>& bin_bounds(std::size_t i) const {
    return bin_sets_[i];
  }

  /**
   * @brief Returns the number of incoming energies / bin boundary sets stored.
   */
  std::size_t size() const { return incoming_energy_.size(); }

 private:
  std::vector<double> incoming_energy_;
  std::vector<std::vector<double>> bin_sets_;

  double sample_bins(double xi1, double xi2,
                     const std::vector<double>& bounds) const;

  double pdf_bins(double E_out, const std::vector<double>& bounds) const;
};

}  // namespace pndl

#endif
