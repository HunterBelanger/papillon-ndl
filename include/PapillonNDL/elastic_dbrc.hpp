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
#ifndef PAPILLON_NDL_ELASTIC_DBRC_H
#define PAPILLON_NDL_ELASTIC_DBRC_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <PapillonNDL/cross_section.hpp>
#include <PapillonNDL/elastic_doppler_broadener.hpp>
#include <functional>
#include <optional>

namespace pndl {

/**
 * @brief This class uses the Doppler Broadening Resonance Correcction
 *        algorithm. It is an improvement on the constant cross section (CXS)
 *        approximation, and provides an exact treatment for the elastic
 *        scattering of neutrons off of nuclides which exhibit strong resonance
 *        behavior at low energies.
 */
class ElasticDBRC : public ElasticDopplerBroadener {
 public:
  /**
   * @param xs The 0 Kelvin elastic scattering cross section for the nuclide.
   */
  ElasticDBRC(const CrossSection& xs) : xs_(xs) {}

  std::array<double, 3> sample_target_velocity(
      const double& Ein, const double& kT, const double& awr,
      const std::function<double()>& rng) const override final;

  std::string algorithm() const override final;

  /**
   * @brief Returns the 0 Kelvin elastic scattering cross section for the
   *        nuclide.
   */
  const CrossSection& elastic_0K_xs() const { return xs_; }

 private:
  CrossSection xs_;

  double max_xs_value(const double& Emin, const double& Emax) const;
};

}  // namespace pndl

#endif
