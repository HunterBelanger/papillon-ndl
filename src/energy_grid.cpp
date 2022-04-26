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
#include <PapillonNDL/energy_grid.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>

namespace pndl {

EnergyGrid::EnergyGrid(const ACE& ace, uint32_t NBINS)
    : energy_values_({0.}), bin_pointers_(), u_min(), du(), urr_start_energy_() {
  energy_values_ = ace.xss(ace.ESZ(), ace.nxs(2));

  // Check if there are URR tables.
  if (ace.jxs(22) != 0) {
    // There are URR tables. Get the starting energy.
    uint32_t URN = static_cast<uint32_t>(ace.jxs(22)) - 1;
    urr_start_energy_ = ace.xss(URN + 6);
  } else {
    // No URR tables. Set start energy to a very large value.
    urr_start_energy_ = 50000;
  }

  if (!std::is_sorted(energy_values_.begin(), energy_values_.end())) {
    std::string mssg =
        "Energy values are not sorted. Index in the XSS block is " +
        std::to_string(ace.ESZ()) + ".";
    throw PNDLException(mssg);
  }

  if (energy_values_.front() <= 0.) {
    std::string mssg =
        "Nevative or zero values in energy grid. Index in the XSS block is " +
        std::to_string(ace.ESZ()) + ".";
    throw PNDLException(mssg);
  }

  hash_energy_grid(NBINS);
}

EnergyGrid::EnergyGrid(const std::vector<double>& energy, uint32_t NBINS)
    : energy_values_(energy), bin_pointers_(), u_min(), du() {
  if (!std::is_sorted(energy_values_.begin(), energy_values_.end())) {
    std::string mssg = "Energy values are not sorted.";
    throw PNDLException(mssg);
  }

  if (energy_values_.front() < 0.) {
    std::string mssg = "Nevative values in energy grid.";
    throw PNDLException(mssg);
  }

  hash_energy_grid(NBINS);
}

void EnergyGrid::hash_energy_grid(uint32_t NBINS) {
  // Generate pointers for lethargy bins
  u_min = std::log(energy_values_.front());
  double u_max = std::log(energy_values_.back());
  du = (u_max - u_min) / static_cast<double>(NBINS);

  bin_pointers_.clear();
  bin_pointers_.reserve(NBINS + 1);

  double E = energy_values_.front();
  std::size_t i = 0;

  // Start by storing index to u_min which is 0
  bin_pointers_.push_back(0);

  // Get energy index for each lethargy bin bound
  for (std::size_t b = 1; b < NBINS + 1; b++) {
    E *= std::exp(du);

    i = std::distance(
            energy_values_.begin(),
            std::lower_bound(energy_values_.begin(), energy_values_.end(), E)) -
        1;

    bin_pointers_.push_back(static_cast<uint32_t>(i));
  }
}

}  // namespace pndl
