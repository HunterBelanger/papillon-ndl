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
#ifndef PAPILLON_NDL_ANGLE_ENERGY_H
#define PAPILLON_NDL_ANGLE_ENERGY_H

/**
 * @file
 * @author Hunter Belanger
 */

#include <functional>
#include <memory>
#include <optional>

namespace pndl {

/**
 * @brief A struct to hold a sampled angle and energy.
 */
struct AngleEnergyPacket {
  double cosine_angle; /**< Sampled cosine of scattering angle */
  double energy;       /**< Sampled outgoing energy in MeV */
};

/**
 * @brief Interface to represent any secondary angle-energy distribution.
 */
class AngleEnergy : public std::enable_shared_from_this<AngleEnergy> {
 public:
  virtual ~AngleEnergy() = default;

  /**
   * @brief Samples an angle and energy from the distribution.
   * @param E_in Incident energy in MeV.
   * @param rng Randum number generation function.
   * @return Sampled cosine of the scattering angle and energy in an
   *         AngleEnergyPacket.
   */
  virtual AngleEnergyPacket sample_angle_energy(
      double E_in, const std::function<double()>& rng) const = 0;

  /**
   * @brief Evaluates the marginal PDF for having a scattering cosine of mu at
   *        incoming energy E_in. Returns an std::optional<double>, as it may
   *        not always be possible to obtain the marginal PDF.
   * @param E_in Incoming energy.
   * @param mu Scattering cosine.
   */
  virtual std::optional<double> angle_pdf(double E_in, double mu) const = 0;

  /**
   * @brief Evaluates the joint PDF for having a scattering cosine of mu at
   *        incoming energy E_in, and exit energy E_out. Returns an
   *        std::optional<double>, as it may not always be possible to
   *        calculate the joint PDF.
   * @param E_in Incoming energy.
   * @param mu Scattering cosine.
   * @param E_out Exit energy.
   */
  virtual std::optional<double> pdf(double E_in, double mu,
                                    double E_out) const = 0;
};

}  // namespace pndl

#endif
