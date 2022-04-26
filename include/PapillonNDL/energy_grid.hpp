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
#ifndef PAPILLON_NDL_ENERGY_GRID_H
#define PAPILLON_NDL_ENERGY_GRID_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/ace.hpp>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

namespace pndl {

/**
 * @brief Holds the hashed energy grid of a Nuclide. An energy grid should
 *        ALWAYS be instantiated as an std::shared_ptr, because a copy is kept
 *        inside all CrossSection instances.
 */
class EnergyGrid : public std::enable_shared_from_this<EnergyGrid> {
  // EnergyGrid inherits from std::enable_shared_from_this, so that
  // Pybind11 can take the references to the energy gird, and form
  // proper shared pointers.

 public:
  /**
   * @param ace ACE file from which to take the energy grid.
   * @param NBINS Number of bins to hash the energy grid into. The
   *              default value is 8192, which is the number of bins
   *              used by MCNP.
   */
  EnergyGrid(const ACE& ace, uint32_t NBINS = 8192);

  /**
   * @param energy Vector of all points in energy grid (sorted).
   * @param NBINS Number of bins to hash the energy grid into. The
   *              default value is 8192, which is the number of bins
   *              used by MCNP.
   */
  EnergyGrid(const std::vector<double>& energy, uint32_t NBINS = 8192);

  ~EnergyGrid() = default;

  /**
   * @brief Returns the ith energy in the grid in MeV.
   * @param i Index into energy grid.
   */
  double operator[](std::size_t i) const { return energy_values_[i]; }

  /**
   * @brief Number of points in the complete energy grid.
   */
  std::size_t size() const { return energy_values_.size(); }

  /**
   * @brief Returns a reference to the energy grid.
   */
  const std::vector<double>& grid() const { return energy_values_; }

  /**
   * @brief Returns the lowest energy in the grid.
   */
  double min_energy() const { return energy_values_.front(); }

  /**
   * @brief Returns the highest energy in the grid.
   */
  double max_energy() const { return energy_values_.back(); }

  /**
   * @brief Returns the starting energy for the unresolved resonance region.
   */
  double urr_min_energy() const { return urr_start_energy_; }

  /**
   * @brief Returns true if the EnergyGrid has an associated unresolved
   *        resonance region.
   */
  bool has_urr() const { return urr_start_energy_ < this->max_energy(); }

  /**
   * @brief Finds the interpolation index for a given energy, using the
   *        hashing algorithm for speed.
   * @param E Energy for which to find the index.
   */
  std::size_t get_lower_index(double E) const {
    if (E <= energy_values_.front()) {
      return 0;
    } else if (E >= energy_values_.back()) {
      return energy_values_.size() - 1;
    }

    // Get current bin
    uint32_t bin = static_cast<uint32_t>((std::log(E) - u_min) / du);

    // lower search index
    uint32_t low_indx = bin_pointers_[bin];
    uint32_t hi_indx = bin_pointers_[bin + 1] + 1;

    std::size_t ind = std::lower_bound(energy_values_.begin() + low_indx,
                                       energy_values_.begin() + hi_indx, E) -
                      energy_values_.begin() - 1;

    return ind;
  }

  /**
   * @brief Re-hashes the energy grid to specified number of pointers.
   * @param NBINS Number of bins to hash the energy grid into.
   */
  void hash_energy_grid(uint32_t NBINS);

 private:
  std::vector<double> energy_values_;
  std::vector<uint32_t> bin_pointers_;
  double u_min, du;
  double urr_start_energy_;
};

}  // namespace pndl

#endif
