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
#ifndef PAPILLON_NDL_ELASTIC_SVT_H
#define PAPILLON_NDL_ELASTIC_SVT_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/elastic_doppler_broadener.hpp>
#include <functional>
#include <memory>
#include <optional>

namespace pndl {

/**
 * @brief This class uses the Sampling of Velocity of Target (SVT) algorithm.
 *        This approach is also sometimes refered to as the Constant Cross
 *        Section (CXS) approximation. It assumes that the microscopic
 *        scattering cross section is approximately constant over the range of
 *        reletive energies which could be observed by the incident neutron (due
 *        to the direction of the target nuclides velcity being random and
 *        isotropically distributed). This method is standard in many Monte
 *        Carlo codes. Despite its ubiquity, it is known to yield inaccurate
 *        results for heavy nuclides which have scattering resonances at low
 *        energies. For these heavy nuclide, the Doppler Broadened Rejection
 *        Correction (DBRC) algorithm is better, as it is an exact treatment.
 */
class ElasticSVT : public ElasticDopplerBroadener {
 public:
  /**
   * @param angle The AngleDistribution for elastic scattering. This
   *              distribution must be given in the center of mass frame.
   * @param awr Atomic weight ratio of the nuclide.
   * @param temperature Temperature in Kelvin of the nuclide.
   * @param use_tar Flag for using the Target At Rest approximation. Default
   *                value is true.
   * @param tar_threshold The threshold for applying the Target At Rest
   *                      approximation. Default value is 400.
   */
  ElasticSVT() = default;

  std::array<double, 3> sample_target_velocity(
      const double& Ein, const double& kT, const double& awr,
      const std::function<double()>& rng) const override final;

  std::string algorithm() const override final;
};

}  // namespace pndl

#endif
