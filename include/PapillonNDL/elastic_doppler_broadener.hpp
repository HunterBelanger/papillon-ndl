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
#ifndef PAPILLON_NDL_ELASTIC_DOPPLER_BROADENER_H
#define PAPILLON_NDL_ELASTIC_DOPPLER_BROADENER_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <array>
#include <functional>
#include <memory>
#include <string>

namespace pndl {

class ElasticDopplerBroadener
    : public std::enable_shared_from_this<ElasticDopplerBroadener> {
 public:
  ElasticDopplerBroadener() = default;
  virtual ~ElasticDopplerBroadener() = default;

  /**
   * @brief Samples the velocity of a target nuclide for using in sampling an
   *        elastic scattering. It is assumed that the target isotope is a
   *        free-gas whose velocity is distributed as a Maxwellian spectrum, at
   *        temperature kT. It is assumed that the direction of the incident
   *        neutron is always along the positive z-axis (0,0,1).
   * @param Ein Incident energy of the neutron in MeV.
   * @param kT Temperature of the "free-gas" in MeV.
   * @param awr Atomic weight ratio of the nuclide.
   * @param rng Random number generator function.
   */
  virtual std::array<double, 3> sample_target_velocity(
      const double& Ein, const double& kT, const double& awr,
      const std::function<double()>& rng) const = 0;

  /**
   * @brief Returns a string with the abbreviation of the elastic kernel
   *        broadening method.
   */
  virtual std::string algorithm() const = 0;
};

}  // namespace pndl

#endif
