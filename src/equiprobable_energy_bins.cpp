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
#include <PapillonNDL/equiprobable_energy_bins.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>
#include <optional>

namespace pndl {

EquiprobableEnergyBins::EquiprobableEnergyBins(const ACE& ace, std::size_t i)
    : incoming_energy_(), bin_sets_() {
  // Get number of interpolation points
  uint32_t NR = ace.xss<uint32_t>(i);
  // Get number of energy points
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);

  // Breakpoints and interpolations are not read, as linear-linear
  // interpolation is always used between incoming energies.

  // Read energies
  incoming_energy_ = ace.xss(i + 2 + 2 * NR, NE);

  if (!std::is_sorted(incoming_energy_.begin(), incoming_energy_.end())) {
    std::string mssg =
        "Incoming energy grid is not sorted. Index of EquiprobableEnergyBins "
        "in XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  uint32_t NET = ace.xss<uint32_t>(i + 2 + 2 * NR + NE);

  // Read energy bins
  for (std::size_t j = 0; j < NE; j++) {
    bin_sets_.push_back(ace.xss(i + 3 + 2 * NR + NE + j * NET, NET));
  }

  // Make sure that each bin set is sorted
  for (std::size_t j = 0; j < bin_sets_.size(); j++) {
    if (!std::is_sorted(bin_sets_[j].begin(), bin_sets_[j].end())) {
      std::string mssg = std::to_string(j) +
                         " bin bounds are not sorted. Index of "
                         "EquiprobableEnergyBins in XSS block is " +
                         std::to_string(i) + ".";
      throw PNDLException(mssg);
    }
  }
}

EquiprobableEnergyBins::EquiprobableEnergyBins(
    const std::vector<double>& incoming_energy,
    const std::vector<std::vector<double>>& bin_bounds)
    : incoming_energy_(incoming_energy), bin_sets_(bin_bounds) {
  if (!std::is_sorted(incoming_energy_.begin(), incoming_energy_.end())) {
    std::string mssg = "Incoming energy grid is not sorted.";
    throw PNDLException(mssg);
  }

  if (incoming_energy_.size() != bin_sets_.size()) {
    std::string mssg =
        "The number of incoming energies does not match the numer of bin sets "
        "provided.";
    throw PNDLException(mssg);
  }

  // Make sure that each bin set is sorted
  for (std::size_t j = 0; j < bin_sets_.size(); j++) {
    if (!std::is_sorted(bin_sets_[j].begin(), bin_sets_[j].end())) {
      std::string mssg = std::to_string(j) + " bin bounds are not sorted.";
      throw PNDLException(mssg);
    }
  }
}

double EquiprobableEnergyBins::sample_energy(
    double E_in, std::function<double()> rng) const {
  // Determine the index of the bounding tabulated incoming energies
  auto in_E_it =
      std::lower_bound(incoming_energy_.begin(), incoming_energy_.end(), E_in);
  if (in_E_it == incoming_energy_.begin()) {
    return sample_bins(rng(), rng(), bin_sets_.front());
  } else if (in_E_it == incoming_energy_.end()) {
    return sample_bins(rng(), rng(), bin_sets_.back());
  }

  std::size_t l = std::distance(incoming_energy_.begin(), in_E_it);
  l--;

  double f = (E_in - incoming_energy_[l]) /
             (incoming_energy_[l + 1] - incoming_energy_[l]);

  if (rng() > f) {
    return sample_bins(rng(), rng(), bin_sets_[l]);
  } else {
    return sample_bins(rng(), rng(), bin_sets_[l + 1]);
  }
}

std::optional<double> EquiprobableEnergyBins::pdf(double E_in,
                                                  double E_out) const {
  // Determine the index of the bounding tabulated incoming energies
  auto in_E_it =
      std::lower_bound(incoming_energy_.begin(), incoming_energy_.end(), E_in);
  if (in_E_it == incoming_energy_.begin()) {
    return pdf_bins(E_out, bin_sets_.front());
  } else if (in_E_it == incoming_energy_.end()) {
    return pdf_bins(E_out, bin_sets_.back());
  }

  std::size_t l = std::distance(incoming_energy_.begin(), in_E_it);
  l--;

  double f = (E_in - incoming_energy_[l]) /
             (incoming_energy_[l + 1] - incoming_energy_[l]);

  return (1. - f) * pdf_bins(E_out, bin_sets_[l]) +
         f * pdf_bins(E_out, bin_sets_[l + 1]);
}

double EquiprobableEnergyBins::sample_bins(
    double xi1, double xi2, const std::vector<double>& bounds) const {
  std::size_t bin = static_cast<std::size_t>(std::floor(bounds.size() * xi1));
  return (bounds[bin + 1] - bounds[bin]) * xi2 + bounds[bin];
}

double EquiprobableEnergyBins::pdf_bins(
    double E_out, const std::vector<double>& bounds) const {
  if (E_out < bounds.front() || E_out > bounds.back()) return 0;

  std::size_t bin = 0;
  for (std::size_t i = 0; i < bounds.size() - 1; i++) {
    if (bounds[i] <= E_out && bounds[i + 1] >= E_out) {
      bin = i;
      break;
    }
  }

  double nbins = static_cast<double>(bounds.size() - 1);
  double prob_per_bin = 1. / nbins;
  double E_low = bounds[bin];
  double E_hi = bounds[bin + 1];
  return prob_per_bin / (E_hi - E_low);
}

}  // namespace pndl
