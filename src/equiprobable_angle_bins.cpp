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
#include <PapillonNDL/equiprobable_angle_bins.hpp>
#include <PapillonNDL/pndl_exception.hpp>
#include <algorithm>
#include <cmath>

namespace pndl {

EquiprobableAngleBins::EquiprobableAngleBins(const ACE& ace, std::size_t i)
    : bounds_(ace.xss(i, NBOUNDS)) {
  if (!std::is_sorted(bounds_.begin(), bounds_.end())) {
    std::string mssg =
        "Bin bounds are not sorted. Index of EquiprobableAngleBins in XSS "
        "block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (bounds_[0] < -1.) {
    std::string mssg =
        "Lowest bin bound is less than -1. Index of EquiprobableAngleBins in "
        "XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg);
  }

  if (bounds_[NBOUNDS - 1] > 1.) {
    std::string mssg =
        "Highest bin bound is more than 1. Index of EquiprobableAngleBins in "
        "XSS block is " +
        std::to_string(i) + ".";
    throw PNDLException(mssg);
  }
}

EquiprobableAngleBins::EquiprobableAngleBins(const std::vector<double>& bounds)
    : bounds_(bounds) {
  if (bounds_.size() != 33) {
    std::string mssg = "Must provide 33 bin boundaries.";
    throw PNDLException(mssg);
  }

  if (!std::is_sorted(bounds_.begin(), bounds_.end())) {
    std::string mssg = "Bin bounds are not sorted.";
    throw PNDLException(mssg);
  }

  if (bounds_[0] < -1.) {
    std::string mssg = "Lowest bin bound is less than -1.";
    throw PNDLException(mssg);
  }

  if (bounds_[NBOUNDS - 1] > 1.) {
    std::string mssg = "Highest bin bound is more than 1.";
    throw PNDLException(mssg);
  }
}

double EquiprobableAngleBins::sample_mu(const std::function<double()>& rng) const {
  const double xi = rng();
  std::size_t bin =
      static_cast<std::size_t>(std::floor(static_cast<double>(NBOUNDS) * xi));
  if (bin == NBOUNDS) bin--;
  double C_b = bin * P_BIN;
  double mu_low = bounds_[bin];
  double mu = ((xi - C_b) / P_BIN) + mu_low;

  if (std::abs(mu) > 1.) mu = std::copysign(1, mu);

  return mu;
}

double EquiprobableAngleBins::pdf(double mu) const {
  if (mu < bounds_.front() || mu > bounds_.back()) return 0.;

  std::size_t bin = 0;
  for (std::size_t i = 0; i < NBOUNDS - 1; i++) {
    if (bounds_[i] <= mu && bounds_[i + 1] >= mu) {
      bin = i;
      break;
    }
  }

  double mu_low = bounds_[bin];
  double mu_hi = bounds_[bin + 1];
  return P_BIN / (mu_hi - mu_low);
}

}  // namespace pndl
