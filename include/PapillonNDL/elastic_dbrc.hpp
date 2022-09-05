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
#ifndef PAPILLON_NDL_ELASTIC_DBRC_H
#define PAPILLON_NDL_ELASTIC_DBRC_H

#include <PapillonNDL/angle_distribution.hpp>
#include <PapillonNDL/angle_energy.hpp>
#include <PapillonNDL/cross_section.hpp>
#include <functional>
#include <optional>

/**
 * @file
 * @author Hunter Belanger
 */

namespace pndl {

/**
 * @brief This class is used to sample elastic scattering of neutrons off of a
 *        nuclide. It uses the Doppler Broadening Resonance Correcction
 *        algorithm. It is an improvement on the constant cross section (CXS)
 *        approximation, and provided an exact treatment for the elastic
 *        scattering of neutrons off of nuclides which exhibit strong resonance
 *        behavior at low energies.
 *
 *        At certain energies, it becomes reasonable to make the approximation
 *        that the target nuclide is at rest, and has no thermal motion. The
 *        threshold for applying this approximation is set with the
 *        tar_threshold parameter. If the incident energy of the neutron (Ein)
 *        is larger than tar_threshold * temperature * k (where k is the
 *        Boltzmann constant), then the target is taken to be stationary. One
 *        exception to this rule is for nuclides with an AWR < 1 (only H1).
 *        Since H1 is actually has a slightly smaller mass than a neutron, the
 *        target at rest approximation is generally inadequate.
 */
class ElasticDBRC : public AngleEnergy {
 public:
  /**
   * @param xs The 0 Kelvin elastic scattering cross section for the nuclide.
   * @param angle The AngleDistribution for elastic scattering. This
   *              distribution must be given in the center of mass frame.
   * @param awr Atomic weight ratio of the nuclide.
   * @param temperature Temperature in Kelvin of the nuclide.
   * @param use_tar Flag for using the Target At Rest approximation. Default
   *                value is true.
   * @param tar_threshold The threshold for applying the Target At Rest
   *                      approximation. Default value is 400.
   */
  ElasticDBRC(const CrossSection& xs, const AngleDistribution& angle,
              double awr, double temperature, bool use_tar = true,
              double tar_threshold = 400.);

  AngleEnergyPacket sample_angle_energy(
      double E_in, std::function<double()> rng) const override final;

  std::optional<double> angle_pdf(double /*E_in*/,
                                  double /*mu*/) const override final {
    return std::nullopt;
  }

  std::optional<double> pdf(double /*E_in*/, double /*mu*/,
                            double /*E_out*/) const override final {
    return std::nullopt;
  }

  /**
   * @brief Returns the AngleDistribution which describes the distribution for
   *        the cosine of the scattering angle in the center-of-mass frame.
   */
  const AngleDistribution& angle_distribution() const { return angle_; }

  /**
   * @brief Returns the 0 Kelvin elastic scattering cross section for the
   *        nuclide.
   */
  const CrossSection& elastic_0K_xs() const { return xs_; }

  /**
   * @breif Returns the Atomic Weight Ratio for the nuclide.
   */
  double awr() const { return awr_; }

  /**
   * @breif Returns the temperature for the nuclide in Kelvin.
   */
  double temperature() const;

  /**
   * @brief If true, the Target At Rest approximation is used for incident
   *        energies which are larger than tar_threshold * kT. If false, the
   *        Target At Rest approximation is never used.
   */
  bool use_tar() const { return use_tar_; };

  /**
   * @brief Returns the threshold for the application of the Target At Rest
   *        approximation.
   */
  double tar_threshold() const { return tar_threshold_; }

 private:
  CrossSection xs_;
  AngleDistribution angle_;
  double awr_;
  double kT_;  // Temperature in MeV
  bool use_tar_;
  double tar_threshold_;
};

}  // namespace pndl

#endif
