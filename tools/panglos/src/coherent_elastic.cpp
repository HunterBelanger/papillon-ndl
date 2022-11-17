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

/**
 * @file
 * @author Hunter Belanger
 */

#include "coherent_elastic.hpp"
#include "interpolator.hpp"

#include <algorithm>
#include <range/v3/to_container.hpp>  // Allows us to use ranges::to<cont<type>>(range);

CoherentElastic::CoherentElastic(const section::Type<7, 2>::CoherentElastic& ce)
    : bragg_edges_(),
      structure_factor_sums_(),
      temperatures_(),
      temp_interps_() {
  bragg_edges_ = ranges::to<std::vector<double>>(ce.energies());
  temperatures_ = ranges::to<std::vector<double>>(ce.temperatures());

  // First, get all the temperature interps
  temp_interps_.reserve(temperatures_.size());
  for (const auto& temp_intrp : ce.temperatureInterpolants()) {
    temp_interps_.push_back(
        Interpolator(static_cast<Interpolation>(temp_intrp)));
  }

  // Get all of the structure factor sums for each temperature
  auto scaterFuncs = ce.S();
  for (std::size_t i = 0; i < temperatures_.size(); i++) {
    structure_factor_sums_.push_back(
        ranges::to<std::vector<double>>(scaterFuncs[i]));
  }

  // Write Information
  std::cout << " Coherent Elastic\n";
  std::cout << " ----------------\n";
  std::cout << " Num. of Bragg Edges = " << bragg_edges_.size() << "\n";
  if (temperature_dependent()) {
    std::cout << " Num. of Temperatures = " << temperatures_.size() << "\n";
  } else {
    std::cout << " No Temperature Dependence\n";
  }
  std::cout << "\n";
}

double CoherentElastic::xs(double T, double Ein) const {
  // First, check if we are bellow the energy threshold
  if (Ein <= bragg_edges_.front()) return 0.;

  // Now, we know we are greater than threshold. Lets find energy index such
  // that bragg_edges_[ie] < Ein < bragg_edges_[ie+1].
  auto eit = std::lower_bound(bragg_edges_.begin(), bragg_edges_.end(), Ein);
  const std::size_t ie = std::distance(bragg_edges_.begin(), eit) - 1;

  // Now we need the S(E,T) value.
  double S = 0.;
  if (temperature_dependent()) {
    // We must interpolate to the correct temperature.
    // First, get the temperature index.
    auto tit = std::lower_bound(temperatures_.begin(), temperatures_.end(), T);

    if (tit == temperatures_.begin()) {
      // T is bellow tabulated temps. We don't extrapolate, just use the
      // lowest provided temperature.
      S = structure_factor_sums_.front()[ie];
    } else if (tit == temperatures_.end()) {
      // T is above the tabulated temps. We don't extrapolate, just use the
      // highest provided temperature.
      S = structure_factor_sums_.back()[ie];
    } else {
      // Interpolate the S value.
      std::size_t it = std::distance(temperatures_.begin(), tit) - 1;
      const double T_low = temperatures_[it];
      const double T_hi = temperatures_[it + 1];
      const double S_low = structure_factor_sums_[it][ie];
      const double S_hi = structure_factor_sums_[it + 1][ie];

      // Interpolate to our temeprature.
      // I think this is the right index. Should probably check.
      S = temp_interps_[it].interpolate(T, T_low, S_low, T_hi, S_hi);
    }
  } else {
    // No temperature dependence. Just use the first set of values.
    S = structure_factor_sums_[0][ie];
  }

  return S / Ein;
}

std::vector<double> CoherentElastic::interpolate_structure_factors(
    double T) const {
  // No temperature dependence, we just return what we have
  if (temperature_dependent() == false) {
    return structure_factor_sums_.front();
  }

  // Desired temperature is outside of the tabulated grid.
  if (T < temperatures_.front() || T > temperatures_.back()) {
  }

  // We must interpolate to the correct temperature.
  // First, get the temperature index.
  auto tit = std::lower_bound(temperatures_.begin(), temperatures_.end(), T);
  if (tit == temperatures_.begin()) {
    // We must be equal to the first temp. Just return that.
    return structure_factor_sums_.front();
  }

  std::vector<double> interpd_struct_fctrs(
      structure_factor_sums_.front().size(), 0.);

  // Interpolate the S value.
  std::size_t it = std::distance(temperatures_.begin(), tit) - 1;
  const double T_low = temperatures_[it];
  const double T_hi = temperatures_[it + 1];

  for (std::size_t i = 0; i < interpd_struct_fctrs.size(); i++) {
    interpd_struct_fctrs[i] =
        temp_interps_[it].interpolate(T, T_low, structure_factor_sums_[it][i],
                                      T_hi, structure_factor_sums_[it + 1][i]);
  }

  return interpd_struct_fctrs;
}
