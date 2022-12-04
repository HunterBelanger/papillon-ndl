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
#include <PapillonNDL/general_evaporation.hpp>
#include <PapillonNDL/tabulated_1d.hpp>
#include <cmath>

namespace pndl {

GeneralEvaporation::GeneralEvaporation(const ACE& ace, std::size_t i)
    : temperature_(), bin_bounds_() {
  uint32_t NR = ace.xss<uint32_t>(i);
  uint32_t NE = ace.xss<uint32_t>(i + 1 + 2 * NR);
  std::vector<uint32_t> NBT;
  std::vector<Interpolation> INT;

  if (NR == 0) {
    NBT = {NE};
    INT = {Interpolation::LinLin};
  } else {
    NBT = ace.xss<uint32_t>(i + 1, NR);
    INT = ace.xss<Interpolation>(i + 1 + NR, NR);
  }

  // Get energy grid
  std::vector<double> energy = ace.xss(i + 2 + 2 * NR, NE);

  std::vector<double> temperature = ace.xss(i + 2 + 2 * NR + NE, NE);

  // Get number of bins
  uint32_t NX = ace.xss<uint32_t>(i + 2 + 2 * NR + 2 * NE);

  // Get bins
  bin_bounds_ = ace.xss(i + 2 + 2 * NR + 2 * NE, NX);

  if (!std::is_sorted(bin_bounds_.begin(), bin_bounds_.end())) {
    std::string mssg =
        "Bin bounds for X are not sorted. Index in the XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  // Create Function1D pointer
  try {
    temperature_ = std::make_shared<Tabulated1D>(NBT, INT, energy, temperature);
  } catch (PNDLException& error) {
    std::string mssg =
        "Could not construct Tabulated1D for the effective nuclear "
        "temperature. "
        "Index in the XSS block is i = " +
        std::to_string(i) + ".";
    error.add_to_exception(mssg);
    throw error;
  }
}

GeneralEvaporation::GeneralEvaporation(std::shared_ptr<Tabulated1D> temperature,
                                       const std::vector<double>& bounds)
    : temperature_(temperature), bin_bounds_(bounds) {
  if (!std::is_sorted(bin_bounds_.begin(), bin_bounds_.end())) {
    std::string mssg = "Bin bounds for X are not sorted.";
    throw PNDLException(mssg);
  }
}

double GeneralEvaporation::sample_energy(
    double E_in, const std::function<double()>& rng) const {
  double T = (*temperature_)(E_in);
  double xi1 = rng();
  std::size_t bin = static_cast<std::size_t>(
      std::floor(static_cast<double>(bin_bounds_.size()) * xi1));
  double xi2 = rng();
  double Chi =
      (bin_bounds_[bin + 1] - bin_bounds_[bin]) * xi2 + bin_bounds_[bin];
  return Chi * T;
}

std::optional<double> GeneralEvaporation::pdf(double E_in, double E_out) const {
  double T = (*temperature_)(E_in);
  double Chi = E_out / T;

  // Go find Chi in bins
  if (Chi < bin_bounds_.front() || Chi > bin_bounds_.back()) return 0.;

  std::size_t bin = 0;
  for (std::size_t i = 0; i < bin_bounds_.size() - 1; i++) {
    if (bin_bounds_[i] <= Chi && bin_bounds_[i + 1] >= Chi) {
      bin = i;
      break;
    }
  }

  double Chi_low = bin_bounds_[bin];
  double Chi_hi = bin_bounds_[bin + 1];
  double nbins = static_cast<double>(bin_bounds_.size() - 1);
  double prob_per_bin = 1. / nbins;

  return prob_per_bin / ((Chi_hi - Chi_low) * T);
}

}  // namespace pndl
